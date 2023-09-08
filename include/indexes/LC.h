//
// Created by joaoleite on 9/5/23.
//

#ifndef GERVLIB_LC_H
#define GERVLIB_LC_H

#include "Index.h"
#include "IndexFactory.h"
#include "NAryTree.h"

namespace gervLib::index::lc
{

    template <typename O, typename T>
    class Node : serialize::Serialize {

    protected:
        size_t nodeID{};
        MEMORY_STATUS memoryStatus{};

    public:
        Node() = default;

        virtual ~Node() = default;

        void setNodeID(size_t _nodeID)
        {
            this->nodeID = _nodeID;
        }

        size_t getNodeID()
        {
            return this->nodeID;
        }

        void setMemoryStatus(MEMORY_STATUS _memoryStatus)
        {
            this->memoryStatus = _memoryStatus;
        }

        MEMORY_STATUS getMemoryStatus()
        {
            return this->memoryStatus;
        }

        virtual bool isLeafNode()
        {
            return false;
        }

        virtual void print(std::ostream& os) const = 0;

        virtual void clear() = 0;

        virtual bool isEqual(std::unique_ptr<Node<O, T>> &other) = 0;

        bool operator==(std::unique_ptr<Node<O, T>> &other)
        {
            return this->isEqual(other);
        }

        bool operator!=(std::unique_ptr<Node<O, T>> &other)
        {
            return !this->isEqual(other);
        }

        std::unique_ptr<u_char[]> serialize() override = 0;

        void deserialize(std::unique_ptr<u_char[]> _data) override = 0;

        size_t getSerializedSize() override = 0;

    };

    template <typename O, typename T>
    class DirectoryNode : public Node<O, T> {

    private:
        std::unique_ptr<Node<O, T>> left, right;
        std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot;
        double r1{}, r2{};

    public:
        DirectoryNode()
        {
            this->left = nullptr;
            this->right = nullptr;
            this->pivot = nullptr;
            this->r1 = std::numeric_limits<double>::max();
            this->r2 = std::numeric_limits<double>::max();
        }

        ~DirectoryNode() override
        {
            if (pivot != nullptr)
            {
                pivot->clear();
                pivot.reset();
            }
        }

        void setLeft(std::unique_ptr<Node<O, T>> _left)
        {
            if (left != nullptr)
                left.reset();

            this->left = std::move(_left);
        }

        void setRight(std::unique_ptr<Node<O, T>> _right)
        {
            if (right != nullptr)
                right.reset();

            this->right = std::move(_right);
        }

        void setPivot(std::unique_ptr<dataset::BasicArrayObject<O, T>> _pivot)
        {
            if (pivot != nullptr)
            {
                pivot->clear();
                pivot.reset();
            }

            this->pivot = std::move(_pivot);
        }

        void setR1(double _range)
        {
            this->r1 = _range;
        }

        void setR2(double _range)
        {
            this->r2 = _range;
        }

        std::unique_ptr<Node<O, T>>& getLeft()
        {
            return this->left;
        }

        std::unique_ptr<Node<O, T>>& getRight()
        {
            return this->right;
        }

        std::unique_ptr<dataset::BasicArrayObject<O, T>>& getPivot()
        {
            return this->pivot;
        }

        double getR1()
        {
            return this->r1;
        }

        double getR2()
        {
            return this->r2;
        }

        bool isLeafNode() override
        {
            return false;
        }

        void print(std::ostream& os) const override
        {
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Memory status: " << this->memoryStatus << std::endl;

            if (pivot != nullptr)
                os << "Pivot: " << *this->pivot << std::endl;
            else
                os << "Pivot: Null" << std::endl;

            os << "R1: " << this->r1 << std::endl;
            os << "R2: " << this->r2 << std::endl;
        }

        void clear()
        {
            if (pivot != nullptr)
            {
                pivot->clear();
                pivot.reset();
            }
        }

        bool isEqual(std::unique_ptr<Node<O, T>> &other)
        {
            if (this->memoryStatus == index::MEMORY_STATUS::IN_DISK)
            {
                return this->nodeID == other->getNodeID();
            }
            else {
                auto *otherDirectoryNode = dynamic_cast<DirectoryNode<O, T>*>(other.get());

                if (otherDirectoryNode == nullptr)
                    return false;

                return (this->pivot->isEqual(*otherDirectoryNode->pivot)) && (this->r1 == otherDirectoryNode->r1) && (this->r2 == otherDirectoryNode->r2);
            }
        }

        std::unique_ptr<u_char[]> serialize()
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(this->getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &this->nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux = index::memoryStatusMap[this->memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;

            if (pivot != nullptr)
            {
                sz = pivot->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> pivotData = pivot->serialize();
                memcpy(data.get() + offset, pivotData.get(), sz);
                offset += sz;
                pivotData.reset();
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            memcpy(data.get() + offset, &this->r1, sizeof(double));
            offset += sizeof(double);

            memcpy(data.get() + offset, &this->r2, sizeof(double));
            offset += sizeof(double);

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data)
        {

            size_t offset = 0, sz;

            memcpy(&this->nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);
            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            this->memoryStatus = index::memoryStatusMapReverse[aux];

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                std::unique_ptr<u_char[]> pivotData = std::make_unique<u_char[]>(sz);
                memcpy(pivotData.get(), _data.get() + offset, sz);
                offset += sz;

                if (pivot != nullptr)
                {
                    pivot->clear();
                    pivot.reset();
                }

                this->pivot = std::make_unique<dataset::BasicArrayObject<O, T>>();
                this->pivot->deserialize(std::move(pivotData));
            }
            else
                this->pivot = nullptr;

            memcpy(&this->r1, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&this->r2, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            _data.reset();

        }

        size_t getSerializedSize()
        {
            size_t ans = 0;

            ans += sizeof(size_t);
            ans += sizeof(size_t) + index::memoryStatusMap[this->memoryStatus].size();
            ans += sizeof(size_t) + (pivot != nullptr ? pivot->getSerializedSize() : 0);
            ans += sizeof(double) * 2;

            return ans;
        }

    };

    template <typename O, typename T>
    class LeafNode : public Node<O, T> {

    private:
        std::unique_ptr<dataset::Dataset<O, T>> dataset;
        std::unique_ptr<Index<O, T>> index;

    public:
        LeafNode()
        {
            this->dataset = nullptr;
            this->index = nullptr;
        }

        ~LeafNode() override
        {
            if (dataset != nullptr)
            {
                dataset->clear();
                dataset.reset();
            }

            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }
        }

        void setDataset(std::unique_ptr<dataset::Dataset<O, T>> _dataset)
        {
            if (dataset != nullptr)
            {
                dataset->clear();
                dataset.reset();
            }

            this->dataset = std::move(_dataset);
        }

        void setIndex(std::unique_ptr<Index<O, T>> _index)
        {
            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }

            this->index = std::move(_index);
        }

        std::unique_ptr<dataset::Dataset<O, T>>& getDataset()
        {
            return this->dataset;
        }

        std::unique_ptr<Index<O, T>>& getIndex()
        {
            return this->index;
        }

        bool isLeafNode() override
        {
            return true;
        }

        void print(std::ostream& os) const override
        {
            os << "Node type: Leaf Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Memory status: " << this->memoryStatus << std::endl;

            if (this->dataset != nullptr)
                os << "Dataset: " << *this->dataset << std::endl;
            else
                os << "Dataset: Null" << std::endl;

            if (this->index != nullptr)
                os << "Index: " << *this->index << std::endl;
            else
                os << "Index: Null" << std::endl;

        }

        void clear() override
        {
            if (dataset != nullptr)
            {
                dataset->clear();
                dataset.reset();
            }

            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }
        }

        bool isEqual(std::unique_ptr<Node<O, T>> &other)
        {
            if (this->memoryStatus == index::MEMORY_STATUS::IN_DISK)
            {
                return this->nodeID == other->getNodeID();
            }
            else {
                auto *otherLeafNode = dynamic_cast<LeafNode<O, T>*>(other.get());

                if (otherLeafNode == nullptr)
                    return false;

                return (this->dataset->isEqual(*otherLeafNode->dataset)) && (this->index->isEqual(otherLeafNode->index));
            }
        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(this->getSerializedSize());
            size_t offset = 0, sz;


            memcpy(data.get() + offset, &this->nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux = index::memoryStatusMap[this->memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;

            if (dataset != nullptr)
            {
                sz = dataset->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> datasetData = dataset->serialize();
                memcpy(data.get() + offset, datasetData.get(), sz);
                offset += sz;
                datasetData.reset();
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (index != nullptr)
            {
                sz = index->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::string indexType = index::indexTypeMap[index->getIndexType()];
                size_t sz2 = indexType.size();

                memcpy(data.get() + offset, &sz2, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, indexType.c_str(), sz2);
                offset += sz2;

                std::unique_ptr<u_char[]> indexData = index->serialize();
                memcpy(data.get() + offset, indexData.get(), sz);
                offset += sz;
                indexData.reset();
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {
            size_t offset = 0, sz;

            memcpy(&this->nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            this->memoryStatus = index::memoryStatusMapReverse[aux];

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                std::unique_ptr<u_char[]> datasetData = std::make_unique<u_char[]>(sz);
                memcpy(datasetData.get(), _data.get() + offset, sz);
                offset += sz;

                if (dataset != nullptr)
                {
                    dataset->clear();
                    dataset.reset();
                }

                this->dataset = std::make_unique<dataset::Dataset<O, T>>();
                this->dataset->deserialize(std::move(datasetData));
            }
            else
                this->dataset = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                size_t sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::string indexType;
                indexType.resize(sz2);

                memcpy(&indexType[0], _data.get() + offset, sz2);
                offset += sz2;

                if (index != nullptr)
                {
                    index->clear();
                    index.reset();
                }

                index = index::IndexFactory<O, T>::createIndex(index::indexTypeMapReverse[indexType]);
                std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
                memcpy(indexData.get(), _data.get() + offset, sz);
                offset += sz;
                index->deserialize(std::move(indexData));
            }
            else
                this->index = nullptr;

            _data.reset();
        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            ans += sizeof(size_t);
            ans += sizeof(size_t) + index::memoryStatusMap[this->memoryStatus].size();
            ans += sizeof(size_t) + (dataset != nullptr ? dataset->getSerializedSize() : 0);

            if (index != nullptr)
            {
                ans += sizeof(size_t) * 2;
                ans += index->getSerializedSize();
                ans += index::indexTypeMap[index->getIndexType()].size();
            }
            else
                ans += sizeof(size_t);

            return ans;
        }

    };

    template <typename O, typename T>
    std::ostream& operator<<(std::ostream& os, const Node<O, T>& printable) {
        printable.print(os);
        return os;
    }

    template <typename O, typename T>
    class LC : public index::Index<O, T> {

    private:
        std::unique_ptr<Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{};
        bool storePivotsInLeaf{}, storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::string serializedTree{};

    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount,iowrite,ioread";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning,iowrite,ioread";
        }

    private:
        enum PositionRelativeToPartition{INSIDE, RING, OUTSIDE};

        PositionRelativeToPartition checkSqPosition(double dist, double mu, double max) {
            if (dist < mu)
                return PositionRelativeToPartition::INSIDE;
            else if (dist <= max)
                return PositionRelativeToPartition::RING;
            else
                return PositionRelativeToPartition::OUTSIDE;
        }

        std::pair<query::Partition<Node<O, T>*>, query::Partition<Node<O, T>*>> getMinMaxDist(double dist, query::Partition<Node<O, T>*> partition)
        {

            auto* directoryNode = (DirectoryNode<O, T>*) partition.getElement();

            PositionRelativeToPartition sqPosition = checkSqPosition(dist, directoryNode->getR1(), directoryNode->getR2());
            query::Partition<Node<O, T>*> left, right;
            left.setElement(directoryNode->getLeft().get());
            right.setElement(directoryNode->getRight().get());

            if (sqPosition == PositionRelativeToPartition::INSIDE)
            {

                //Min left
                left.setMin(0.0);

                //Max left
                if (partition.getMax() > directoryNode->getR1() + dist)
                    left.setMax(directoryNode->getR1() + dist);
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(directoryNode->getR1() - dist);

                //Max right
                if (partition.getMax() > directoryNode->getR2() + dist)
                    right.setMax(directoryNode->getR2() + dist);
                else
                    right.setMax(partition.getMax());

            }
            else if (sqPosition == PositionRelativeToPartition::RING)
            {

                //Min left
                left.setMin(dist - directoryNode->getR1());

                //Max left
                if (partition.getMax() > dist + directoryNode->getR1())
                    left.setMax(dist + directoryNode->getR1());
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(0.0);

                //Max right
                if (partition.getMax() > directoryNode->getR2() + dist)
                    right.setMax(directoryNode->getR2() + dist);
                else
                    right.setMax(partition.getMax());

            }
            else
            {

                //Min left
                left.setMin(dist - directoryNode->getR1());

                //Max left
                if (partition.getMax() > dist + directoryNode->getR1())
                    left.setMax(dist + directoryNode->getR1());
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(dist - directoryNode->getR2());

                //Max right
                if (partition.getMax() > dist + directoryNode->getR2())
                    right.setMax(dist + directoryNode->getR2());
                else
                    right.setMax(partition.getMax());

            }

            return std::make_pair(left, right);

        }

        void clearRecursive(std::unique_ptr<Node<O, T>>& node)
        {
            if (node == nullptr)
                return;

            if (!node->isLeafNode())
            {
                auto* directoryNode = (DirectoryNode<O, T>*) node.get();
                clearRecursive(directoryNode->getLeft());
                clearRecursive(directoryNode->getRight());
            }

            node->clear();

        }

        void deleteRecursive(std::unique_ptr<Node<O, T>> node)
        {
            if (node == nullptr)
                return;

            if (!node->isLeafNode())
            {
                auto* directoryNode = (DirectoryNode<O, T>*) node.get();
                deleteRecursive(std::move(directoryNode->getLeft()));
                deleteRecursive(std::move(directoryNode->getRight()));
            }

            node->clear();
            node.reset();

        }

        bool isEqualHelper(std::unique_ptr<Node<O, T>>& node1, std::unique_ptr<Node<O, T>>& node2)
        {
            if (node1 == nullptr && node2 == nullptr)
                return true;

            if ((node1 != nullptr && node2 == nullptr) || (node1 == nullptr && node2 != nullptr))
                return false;

            if (!node1->isEqual(node2))
                return false;

            if ((!node1->isLeafNode() && node2->isLeafNode()) || (node1->isLeafNode() && !node2->isLeafNode()))
                return false;
            else if (!node1->isLeafNode() && !node2->isLeafNode())
            {
                auto* directoryNode1 = (DirectoryNode<O, T>*) node1.get();
                auto* directoryNode2 = (DirectoryNode<O, T>*) node2.get();

                return isEqualHelper(directoryNode1->getLeft(), directoryNode2->getLeft()) && isEqualHelper(directoryNode1->getRight(), directoryNode2->getRight());
            }
            else
                return true;

        }

        std::string serializeTreeRecursive(std::unique_ptr<Node<O, T>>& node)
        {
            if (node == nullptr)
                return "null ";

            if (node->getMemoryStatus() != index::MEMORY_STATUS::IN_DISK)
            {
                std::unique_ptr<u_char[]> data = node->serialize();
                this->pageManager->save(node->getNodeID(), std::move(data), node->getSerializedSize());
                node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);
            }

            std::string result = (node->isLeafNode() ? "L" : "D") + std::to_string(node->getNodeID()) + " 2 ";

            if (node->isLeafNode())
            {
                result += "null ";
                result += "null ";
            }
            else{
                auto* directoryNode = (DirectoryNode<O, T>*) node.get();
                result += serializeTreeRecursive(directoryNode->getLeft());
                result += serializeTreeRecursive(directoryNode->getRight());
            }

            return result;

        }

        std::unique_ptr<naryTree::NodeNAry> deserializeTreeRecursive(std::stringstream& ss)
        {
            std::string valStr;
            ss >> valStr;

            if (valStr == "null")
                return nullptr;

            size_t numChildren;
            ss >> numChildren;

            std::unique_ptr<naryTree::NodeNAry> node = std::make_unique<naryTree::NodeNAry>(valStr);

            for (size_t i = 0; i < numChildren; i++)
                node->children.push_back(deserializeTreeRecursive(ss));

            return node;
        }

        void buildTree(std::unique_ptr<Node<O, T>>& node, std::unique_ptr<naryTree::NodeNAry>& aux)
        {

            if (aux == nullptr)
                return;

            if (aux->value[0] == 'L')
            {
                node = std::make_unique<LeafNode<O, T>>();
                node->setNodeID(std::stoull(aux->value.substr(1)));

                if (!storeLeafNode) {
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_MEMORY);
                    std::unique_ptr<u_char[]> data = this->pageManager->load(node->getNodeID());
                    node->deserialize(std::move(data));
                }
                else
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);

            }
            else
            {
                node = std::make_unique<DirectoryNode<O, T>>();
                node->setNodeID(std::stoull(aux->value.substr(1)));

                if (!storeDirectoryNode) {
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_MEMORY);
                    std::unique_ptr<u_char[]> data = this->pageManager->load(node->getNodeID());
                    node->deserialize(std::move(data));
                }
                else
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);

                auto* directoryNode = (DirectoryNode<O, T>*) node.get();

                buildTree(directoryNode->getLeft(), aux->children[0]);
                buildTree(directoryNode->getRight(), aux->children[1]);

            }

        }

    public:
        LC()
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->root = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->storePivotsInLeaf = false;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = index::INDEX_TYPE::LC;
            this->indexName = "LC";
            this->indexFolder = "";
        }

        LC(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
           std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
           std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0,
           bool _storePivotsInLeaf = true, bool _storeDirectoryNode = false, bool _storeLeafNode = false, bool _useLAESA = true, std::string folder="")
        {
            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->root = nullptr;
            this->pageSize = _pageSize;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->storePivotsInLeaf = _storePivotsInLeaf;
            this->numPerLeaf = _numPerLeaf;
            this->numPivots = _numPivots;
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->useLAESA = _useLAESA;
            this->indexType = index::INDEX_TYPE::LC;
            this->indexName = "LC";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("lc_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit LC(std::string _folder, std::string serializedFile = "")
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->root = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->storePivotsInLeaf = false;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = index::INDEX_TYPE::LC;
            this->indexName = "LC";

            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);
        }

        ~LC()
        {
            deleteRecursive(std::move(root));
        }

        //TODO implement serialize, deserialize, getSerializedSize

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            size_t currentNodeID = 0;
            std::queue<std::pair<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>>> nodeQueue;
            std::pair<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>> currentNode;
            std::unique_ptr<query::Result<dataset::BasicArrayObject<O, T>>> result = std::make_unique<query::Result<dataset::BasicArrayObject<O, T>>>();
            result->setMaxSize(this->numPerLeaf);
            std::unique_ptr<pivots::Pivot<O, T>> globalPivots;
            query::ResultEntry<dataset::BasicArrayObject<O, T>> entry;
            double dist;
            std::vector<std::pair<double, size_t>> distances;

            if (useLAESA)
            {
                globalPivots = pivots::PivotFactory<O, T>::createPivot(this->pivots->getPivotType(), this->pivots);
                globalPivots->operator()(this->dataset, this->distanceFunction, this->numPivots);
            }

            if (this->dataset->getCardinality() <= numPerLeaf)
                root = std::make_unique<LeafNode<O, T>>();
            else
                root = std::make_unique<DirectoryNode<O, T>>();

            nodeQueue.push(std::make_pair(root.get(), std::move(this->dataset)));

            while (!nodeQueue.empty())
            {

                currentNode = std::move(nodeQueue.front());
                nodeQueue.pop();

                if (currentNode.second->getCardinality() <= numPerLeaf)
                {

                    auto* leafNode = (LeafNode<O, T>*) currentNode.first;
                    leafNode->setNodeID(currentNodeID++);

                    if (useLAESA) {

                        std::filesystem::path leafIndexPath(this->indexFolder);
                        leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                this->distanceFunction->getDistanceType());
                        size_t oldDistCount = df->getDistanceCount();
                        std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                std::move(currentNode.second), std::move(df), pivots::PivotFactory<O, T>::clone(globalPivots), this->numPivots,
                                leafIndexPath);
                        idx->getDistanceFunction()->setDistanceCount(idx->getDistanceFunction()->getDistanceCount() + oldDistCount);
                        leafNode->setIndex(std::move(idx));

                    }
                    else
                        leafNode->setDataset(std::move(currentNode.second));

                    if (storeLeafNode)
                    {
                        leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> leafData = leafNode->serialize();
                        this->pageManager->save(leafNode->getNodeID(), std::move(leafData), leafNode->getSerializedSize());
                        leafNode->clear();
                    }
                    else
                        leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                }
                else
                {

                    auto *directoryNode = (DirectoryNode<O, T> *) currentNode.first;
                    directoryNode->setNodeID(currentNodeID++);

                    this->pivots->operator()(currentNode.second, this->distanceFunction, 1);

                    if (!storePivotsInLeaf)
                        currentNode.second->erase(this->pivots->getPivot(0));

                    for (size_t i = 0; i < currentNode.second->getCardinality(); i++)
                    {

                        dist = this->distanceFunction->operator()(currentNode.second->getElement(i), this->pivots->getPivot(0));
                        distances.emplace_back(dist, i);

                    }

                    std::sort(distances.begin(), distances.end(), [](std::pair<double, size_t> a, std::pair<double, size_t> b) {
                        return a.first < b.first;
                    });

                    std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot = std::make_unique<dataset::BasicArrayObject<O, T>>(this->pivots->getPivot(0));
                    directoryNode->setPivot(std::move(pivot));
                    directoryNode->setR1(distances[numPerLeaf-1].first);
                    directoryNode->setR2(distances.back().first);

                    std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>(), rightDataset = std::make_unique<dataset::Dataset<O, T>>();
                    leftDataset->setPath(currentNode.second->getPath());
                    rightDataset->setPath(currentNode.second->getPath());
                    leftDataset->setSeed(currentNode.second->getSeed());
                    rightDataset->setSeed(currentNode.second->getSeed());

                    for (size_t i = 0; i < distances.size(); i++)
                    {
                        if (i < numPerLeaf)
                            leftDataset->insert(currentNode.second->getElement(distances[i].second));
                        else
                            rightDataset->insert(currentNode.second->getElement(distances[i].second));
                    }

                    if (leftDataset->getCardinality() <= numPerLeaf)
                    {
                        std::unique_ptr<lc::Node<O, T>> leftLeafNode = std::make_unique<LeafNode<O, T>>();
                        nodeQueue.push(std::make_pair(leftLeafNode.get(), std::move(leftDataset)));
                        directoryNode->setLeft(std::move(leftLeafNode));
                    }
                    else
                    {
                        std::unique_ptr<lc::Node<O, T>> leftDirectoryNode = std::make_unique<DirectoryNode<O, T>>();
                        nodeQueue.push(std::make_pair(leftDirectoryNode.get(), std::move(leftDataset)));
                        directoryNode->setLeft(std::move(leftDirectoryNode));
                    }

                    if (rightDataset->getCardinality() <= numPerLeaf)
                    {
                        std::unique_ptr<lc::Node<O, T>> rightLeafNode = std::make_unique<LeafNode<O, T>>();
                        nodeQueue.push(std::make_pair(rightLeafNode.get(), std::move(rightDataset)));
                        directoryNode->setRight(std::move(rightLeafNode));
                    }
                    else
                    {
                        std::unique_ptr<lc::Node<O, T>> rightDirectoryNode = std::make_unique<DirectoryNode<O, T>>();
                        nodeQueue.push(std::make_pair(rightDirectoryNode.get(), std::move(rightDataset)));
                        directoryNode->setRight(std::move(rightDirectoryNode));
                    }

                    if (this->storeDirectoryNode) {
                        directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                        this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                        directoryNode->clear();
                    }
                    else
                        directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                    currentNode.second->clear();
                    currentNode.second.reset();


                    distances.clear();

                }

            }

            timer.stop();

            std::ofstream buildFile(this->buildFile, std::ios::app);
            buildFile << timer.getElapsedTime()
                      << ","
                      << timer.getElapsedTimeSystem()
                      << ","
                      << timer.getElapsedTimeUser()
                      << ","
                      << this->distanceFunction->getDistanceCount()
                      << ","
                      << std::to_string(configure::IOWrite - ioW)
                      << ","
                      << std::to_string(configure::IORead - ioR) << std::endl;
            buildFile.close();

        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<LC<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

        }

        std::unique_ptr<Node<O, T>>& getRoot()
        {
            return root;
        }

        void print(std::ostream& os) const override {

            std::stack<Node<O, T>*> nodeStack;
            nodeStack.push(root.get());

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "LC: " << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Store pivots in leaf: " << (storePivotsInLeaf ? "true" : "false") << std::endl;
            os << "Number of objects per leaf: " << numPerLeaf << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;
            os << "Use LAESA: " << (useLAESA ? "true" : "false") << std::endl;

            while (!nodeStack.empty())
            {

                Node<O, T>* currentNode = nodeStack.top();
                nodeStack.pop();

                os << "**********************************************************************************************************************************************************************************\n\n";
                os << *currentNode << std::endl;

                if (!currentNode->isLeafNode())
                {
                    auto *directoryNode = (DirectoryNode<O, T>*) currentNode;
                    nodeStack.push(directoryNode->getLeft().get());
                    nodeStack.push(directoryNode->getRight().get());
                }

            }

        }

        void clear() override
        {
            clearRecursive(root);
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults, bool saveStatistics) override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            this->leafNodeAccess = 0;
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            std::priority_queue<query::Partition<Node<O, T>*>, std::vector<query::Partition<Node<O, T>*>>, std::greater<query::Partition<Node<O, T>*>>> nodeQueue;
            std::priority_queue<query::ResultEntry<O>, std::vector<query::ResultEntry<O>>, std::greater<query::ResultEntry<O>>> elementQueue;
            query::Result<O> result;
            result.setMaxSize(k);
            query::Partition<Node<O, T>*> currentPartition;
            Node<O, T>* currentNode;
            DirectoryNode<O, T>* currentDirectoryNode;
            LeafNode<O, T>* currentLeafNode;
            LAESA<O, T>* laesa;
            double dist;

            nodeQueue.push(query::Partition<Node<O, T>*>(root.get(), 0.0, std::numeric_limits<double>::max()));

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty()))
            {

                if (elementQueue.empty())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if (currentNode->isLeafNode())
                    {

                        currentLeafNode = (LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                    }
                    else
                    {
                        currentDirectoryNode = (DirectoryNode<O, T>*) currentNode;

                        dist = this->distanceFunction->operator()(query, *currentDirectoryNode->getPivot());

                        if (!storePivotsInLeaf)
                            elementQueue.push(query::ResultEntry<O>(currentDirectoryNode->getPivot()->getOID(), dist));

                        std::pair<query::Partition<Node<O, T>*>, query::Partition<Node<O, T>*>> partitions = getMinMaxDist(dist, currentPartition);
                        nodeQueue.push(partitions.first);
                        nodeQueue.push(partitions.second);
                    }

                    if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                    {
                        currentNode->clear();
                    }

                }
                else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if (currentNode->isLeafNode())
                    {

                        currentLeafNode = (LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                    }
                    else
                    {

                        currentDirectoryNode = (DirectoryNode<O, T>*) currentNode;

                        dist = this->distanceFunction->operator()(query, *currentDirectoryNode->getPivot());

                        if (!storePivotsInLeaf)
                            elementQueue.push(query::ResultEntry<O>(currentDirectoryNode->getPivot()->getOID(), dist));

                        std::pair<query::Partition<Node<O, T>*>, query::Partition<Node<O, T>*>> partitions = getMinMaxDist(dist, currentPartition);
                        nodeQueue.push(partitions.first);
                        nodeQueue.push(partitions.second);
                    }

                    if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                    {
                        currentNode->clear();
                    }

                }
                else
                {

                    result.push(elementQueue.top());
                    elementQueue.pop();

                }

            }

            std::vector<query::ResultEntry<O>> ans = result.getResults();
            std::reverse(ans.begin(), ans.end());

            std::string expt_id = utils::generateExperimentID();
            timer.stop();

            if (saveResults)
            {
                this->saveResultToFile(ans, query, "kNNIncremental", expt_id);
            }

            if (saveStatistics)
            {
                this->saveStatistics({expt_id, std::to_string(k), "-1",
                                      std::to_string(timer.getElapsedTime()),
                                      std::to_string(timer.getElapsedTimeSystem()),
                                      std::to_string(timer.getElapsedTimeUser()),
                                      std::to_string(this->distanceFunction->getDistanceCount()),
                                      std::to_string(this->prunning),
                                      std::to_string(configure::IOWrite - ioW),
                                      std::to_string(configure::IORead - ioR)});
            }

            return ans;

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(this->getSerializedSize());
            size_t offset = 0, sz;

            sz = gervLib::index::Index<O, T>::getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = gervLib::index::Index<O, T>::serialize();
            memcpy(data.get() + offset, indexData.get(), sz);
            offset += sz;
            indexData.reset();

            memcpy(data.get() + offset, &numPerLeaf, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numPivots, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storePivotsInLeaf ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storeDirectoryNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storeLeafNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (useLAESA ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = serializedTree.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, serializedTree.c_str(), sz);
            offset += sz;

            serializedTree.clear();

            return data;
        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {
            size_t offset = 0, sz;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
            memcpy(indexData.get(), _data.get() + offset, sz);
            offset += sz;
            gervLib::index::Index<O, T>::deserialize(std::move(indexData));

            memcpy(&numPerLeaf, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numPivots, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storePivotsInLeaf = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeDirectoryNode = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeLeafNode = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            useLAESA = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            std::stringstream ss(aux);
            auto root_helper = deserializeTreeRecursive(ss);

            if (root != nullptr)
                root.reset();

            buildTree(root, root_helper);

            _data.reset();
        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            this->serializedTree = serializeTreeRecursive(root);

            ans += sizeof(size_t) + gervLib::index::Index<O, T>::getSerializedSize();
            ans += sizeof(size_t) * 6;

            ans += sizeof(size_t) + serializedTree.size();

            return ans;
        }

    };

}

#endif //GERVLIB_LC_H
