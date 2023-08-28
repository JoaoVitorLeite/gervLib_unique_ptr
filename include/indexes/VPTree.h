//
// Created by joaovictor on 21/08/23.
//

#ifndef GERVLIB_VPTREE_H
#define GERVLIB_VPTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "PivotFactory.h"
#include "DistanceFactory.h"
#include "NAryTree.h"
#include <queue>

namespace gervLib::index::vptree
{

    template <typename O, typename T>
    class Node : serialize::Serialize
    {

    protected:
        std::unique_ptr<Node<O, T>> left, right;
        std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot;
        double mu{}, coverage{};
        index::MEMORY_STATUS memoryStatus = index::MEMORY_STATUS::NONE;
        size_t nodeID{};

    public:
        Node() = default;

        virtual ~Node()
        {
            if (this->pivot != nullptr)
                this->pivot.reset();
        }

        virtual void clear()
        {
            if (this->pivot != nullptr)
                this->pivot.reset();
        }

        void setPivot(std::unique_ptr<dataset::BasicArrayObject<O, T>> _pivot)
        {
            clear();
            this->pivot = std::move(_pivot);
        }

        void setMu(double _mu)
        {
            this->mu = _mu;
        }

        void setCoverage(double _coverage)
        {
            this->coverage = _coverage;
        }

        void setMemoryStatus(index::MEMORY_STATUS _memoryStatus)
        {
            this->memoryStatus = _memoryStatus;
        }

        void setNodeID(size_t _nodeID)
        {
            this->nodeID = _nodeID;
        }

        std::unique_ptr<dataset::BasicArrayObject<O, T>> &getPivot()
        {
            return this->pivot;
        }

        double getMu()
        {
            return this->mu;
        }

        double getCoverage()
        {
            return this->coverage;
        }

        index::MEMORY_STATUS getMemoryStatus()
        {
            return this->memoryStatus;
        }

        size_t getNodeID()
        {
            return this->nodeID;
        }

        virtual bool isEqual(std::unique_ptr<Node<O, T>> &other)
        {
            if (memoryStatus == index::MEMORY_STATUS::IN_DISK)
            {
                return nodeID == other->nodeID;
            }
            else
                return (this->pivot->isEqual(*other->pivot)) && (this->mu == other->mu) && (this->coverage == other->coverage);
        }

        bool operator==(std::unique_ptr<Node<O, T>> &other)
        {
            return this->isEqual(other);
        }

        bool operator!=(std::unique_ptr<Node<O, T>> &other)
        {
            return !this->isEqual(other);
        }

        virtual bool isLeafNode()
        {
            return false;
        }

        virtual void print(std::ostream& os) const
        {

            os << "Node type = Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;
            if (pivot != nullptr)
                os << "Pivot: " << *pivot << std::endl;
            else
                os << "Pivot: NULL" << std::endl;
            os << "Mu: " << mu << std::endl;
            os << "Coverage: " << coverage << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[memoryStatus] << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            sz = (pivot == nullptr ? 0 : pivot->getSerializedSize());
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            if (pivot != nullptr)
            {
                std::unique_ptr<u_char[]> pivotData = pivot->serialize();
                memcpy(data.get() + offset, pivotData.get(), sz);
                offset += sz;
                pivotData.reset();
            }

            memcpy(data.get() + offset, &mu, sizeof(double));
            offset += sizeof(double);

            memcpy(data.get() + offset, &coverage, sizeof(double));
            offset += sizeof(double);

            memcpy(data.get() + offset, &nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux = index::memoryStatusMap[memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                std::unique_ptr<u_char[]> pivotData = std::make_unique<u_char[]>(sz);
                memcpy(pivotData.get(), _data.get() + offset, sz);
                offset += sz;
                clear();
                pivot = std::make_unique<dataset::BasicArrayObject<O, T>>();
                pivot->deserialize(std::move(pivotData));
            }

            memcpy(&mu, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&coverage, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            memoryStatus = index::memoryStatusMapReverse[aux];

            _data.reset();

        }

        size_t getSerializedSize() override
        {

            return sizeof(size_t) + (pivot == nullptr ? 0 : pivot->getSerializedSize()) + //pivot
                   sizeof(double) * 2 + //mu and coverage
                   sizeof(size_t) + //nodeID
                   sizeof(size_t) + gervLib::index::memoryStatusMap[memoryStatus].size(); //memoryStatus

        }

        std::unique_ptr<Node<O, T>> &getLeft()
        {
            return this->left;
        }

        std::unique_ptr<Node<O, T>> &getRight()
        {
            return this->right;
        }

        void setLeft(std::unique_ptr<Node<O, T>> _left)
        {
            this->left.reset();
            this->left = std::move(_left);
        }

        void setRight(std::unique_ptr<Node<O, T>> _right)
        {
            this->right.reset();
            this->right = std::move(_right);
        }

    };

    template <typename O, typename T>
    class DirectoryNode : public Node<O, T>
    {

    public:
        DirectoryNode(): Node<O, T>()
        {
            this->left = nullptr;
            this->right = nullptr;
        }

        DirectoryNode(std::unique_ptr<Node<O, T>> _left, std::unique_ptr<Node<O, T>> _right): Node<O, T>()
        {
            this->left = std::move(_left);
            this->right = std::move(_right);
        }

        DirectoryNode(std::unique_ptr<dataset::BasicArrayObject<O, T>> _pivot, std::unique_ptr<Node<O, T>> _left, std::unique_ptr<Node<O, T>> _right): Node<O, T>()
        {
            this->setPivot(std::move(_pivot));
            this->left = std::move(_left);
            this->right = std::move(_right);
        }

        ~DirectoryNode() override = default;

        bool isLeafNode() override
        {
            return false;
        }

        void print(std::ostream& os) const override
        {

            os << "Node type = Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            if (this->pivot != nullptr)
                os << "Pivot: " << *this->pivot << std::endl;
            else
                os << "Pivot: NULL" << std::endl;
            os << "Mu: " << this->mu << std::endl;
            os << "Coverage: " << this->coverage << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;

        }

    };

    template <typename O, typename T>
    class LeafNode : public Node<O, T>
    {

    private:
        std::unique_ptr<dataset::Dataset<O, T>> dataset;
        std::unique_ptr<Index<O, T>> index;

    public:
        LeafNode() : Node<O, T>()
        {
            this->dataset = nullptr;
            this->index = nullptr;
        }

        explicit LeafNode(std::unique_ptr<dataset::Dataset<O, T>> _dataset) : Node<O, T>()
        {
            this->dataset = std::move(_dataset);
            this->index = nullptr;
        }


        LeafNode(std::unique_ptr<dataset::Dataset<O, T>> _dataset, std::unique_ptr<Index<O, T>> _index) : Node<O, T>()
        {
            this->dataset = std::move(_dataset);
            this->index = std::move(_index);
        }

        ~LeafNode() override
        {
            Node<O, T>::clear();
            if (this->dataset != nullptr)
                this->dataset.reset();
            if (this->index != nullptr)
                this->index.reset();
        }

        bool isLeafNode() override
        {
            return true;
        }

        void clear() override
        {
            Node<O, T>::clear();
            if (this->dataset != nullptr)
                this->dataset.reset();
            if (this->index != nullptr)
                this->index.reset();
        }

        void setIndex(std::unique_ptr<Index<O, T>> _index)
        {
            this->index.reset();
            this->index = std::move(_index);
        }

        void setDataset(std::unique_ptr<dataset::Dataset<O, T>> _dataset)
        {
            this->dataset.reset();
            this->dataset = std::move(_dataset);
        }

        std::unique_ptr<Index<O, T>> &getIndex()
        {
            return this->index;
        }

        std::unique_ptr<dataset::Dataset<O, T>> &getDataset()
        {
            return this->dataset;
        }

        void insert(std::unique_ptr<dataset::BasicArrayObject<O, T>> &object)
        {
            if (dataset == nullptr)
                dataset = std::make_unique<dataset::Dataset<O, T>>();
            this->dataset->insert(*object);
        }

        bool isEqual(std::unique_ptr<Node<O, T>> &other) override
        {
            auto* node = dynamic_cast<LeafNode<O, T>*>(other.get());

            if (this->memoryStatus == index::MEMORY_STATUS::IN_DISK && node->memoryStatus == index::MEMORY_STATUS::IN_DISK)
                return this->nodeID == node->nodeID;

            if (node == nullptr)
                return false;

            if ((this->dataset == nullptr && node->dataset != nullptr) || (this->dataset != nullptr && node->dataset == nullptr))
                return false;

            if (this->dataset != nullptr && node->dataset != nullptr)
                if (!this->dataset->isEqual(*node->dataset))
                    return false;

            if ((this->index == nullptr && node->index != nullptr) || (this->index != nullptr && node->index == nullptr))
                return false;

            if (this->index != nullptr && node->index != nullptr)
                if (!this->index->isEqual(node->index))
                    return false;

            if ((this->mu != node->mu) && (this->coverage != node->coverage))
                return false;

            if ((this->pivot == nullptr && node->pivot != nullptr) || (this->pivot != nullptr && node->pivot == nullptr))
                return false;

            if (this->pivot != nullptr && node->pivot != nullptr)
                if (!this->pivot->isEqual(*node->pivot))
                    return false;

            return true;

        }

        void print(std::ostream& os) const override
        {

            os << "Node type = Leaf Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            if (this->pivot != nullptr)
                os << "Pivot: " << *this->pivot << std::endl;
            else
                os << "Pivot: NULL" << std::endl;
            os << "Mu: " << this->mu << std::endl;
            os << "Coverage: " << this->coverage << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;

            if (this->dataset != nullptr)
                os << "Dataset: " << *this->dataset << std::endl;
            else
                os << "Dataset: NULL" << std::endl;

            if (this->index != nullptr)
                os << "Index: " << *this->index << std::endl;
            else
                os << "Index: NULL" << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            sz = Node<O, T>::getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> nodeData = Node<O, T>::serialize();
            memcpy(data.get() + offset, nodeData.get(), sz);
            offset += sz;
            nodeData.reset();

            sz = (dataset == nullptr ? 0 : dataset->getSerializedSize());
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            if (dataset != nullptr)
            {
                std::unique_ptr<u_char[]> datasetData = dataset->serialize();
                memcpy(data.get() + offset, datasetData.get(), sz);
                offset += sz;
                datasetData.reset();
            }

            sz = (index == nullptr ? 0 : index->getSerializedSize());
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            if (index != nullptr)
            {
                std::string aux = index::indexTypeMap[index->getIndexType()];
                size_t sz2 = aux.size();

                memcpy(data.get() + offset, &sz2, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, aux.c_str(), sz2);
                offset += sz2;

                std::unique_ptr<u_char[]> indexData = index->serialize();
                memcpy(data.get() + offset, indexData.get(), sz);
                offset += sz;
                indexData.reset();

            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> nodeData = std::make_unique<u_char[]>(sz);
            memcpy(nodeData.get(), _data.get() + offset, sz);
            offset += sz;
            Node<O, T>::deserialize(std::move(nodeData));

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                std::unique_ptr<u_char[]> datasetData = std::make_unique<u_char[]>(sz);
                memcpy(datasetData.get(), _data.get() + offset, sz);
                offset += sz;
                dataset.reset();
                dataset = std::make_unique<dataset::Dataset<O, T>>();
                dataset->deserialize(std::move(datasetData));
            }

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                size_t sz2;

                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::string aux;
                aux.resize(sz2);

                memcpy(&aux[0], _data.get() + offset, sz2);
                offset += sz2;

                index.reset();
                index = index::IndexFactory<O, T>::createIndex(index::indexTypeMapReverse[aux]);

                std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
                memcpy(indexData.get(), _data.get() + offset, sz);
                offset += sz;
                index->deserialize(std::move(indexData));

            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = sizeof(size_t) + Node<O, T>::getSerializedSize() + // node
                         sizeof(size_t) + (dataset != nullptr ? dataset->getSerializedSize() : 0); // Dataset

            if (index != nullptr)
            {
                ans += sizeof(size_t) * 2 + index::indexTypeMap[index->getIndexType()].size() + index->getSerializedSize();
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
    class VPTree : public index::Index<O, T>
    {

    private:
        std::unique_ptr<Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{};
        bool storePivotsInLeaf{}, storeLeafNode{}, storeDirectoryNode{};
        std::string serializedTree{};

    private:
        void deleteRecursive(std::unique_ptr<Node<O, T>> node)
        {
            if (node == nullptr)
                return;

            deleteRecursive(std::move(node->getLeft()));
            deleteRecursive(std::move(node->getRight()));
            node->clear();
            node.reset();

        }

        void clearRecursive(std::unique_ptr<Node<O, T>>& node)
        {
            if (node == nullptr)
                return;

            clearRecursive(node->getLeft());
            clearRecursive(node->getRight());
            node->clear();

        }

        bool isEqualHelper(std::unique_ptr<Node<O, T>>& node1, std::unique_ptr<Node<O, T>>& node2)
        {
            if (node1 == nullptr && node2 == nullptr)
                return true;

            if ((node1 != nullptr && node2 == nullptr) || (node1 == nullptr && node2 != nullptr))
                return false;

            if (!node1->isEqual(node2))
                return false;

            return isEqualHelper(node1->getLeft(), node2->getLeft()) && isEqualHelper(node1->getRight(), node2->getRight());
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
            result += serializeTreeRecursive(node->getLeft());
            result += serializeTreeRecursive(node->getRight());

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

                node->setLeft(nullptr);
                node->setRight(nullptr);
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

                node->setLeft(nullptr);
                node->setRight(nullptr);
            }

            buildTree(node->getLeft(), aux->children[0]);
            buildTree(node->getRight(), aux->children[1]);

        }

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
            PositionRelativeToPartition sqPosition = checkSqPosition(dist, partition.getElement()->getMu(), partition.getElement()->getCoverage());
            query::Partition<Node<O, T>*> left, right;
            left.setElement(partition.getElement()->getLeft().get());
            right.setElement(partition.getElement()->getRight().get());

            if (sqPosition == PositionRelativeToPartition::INSIDE)
            {

                //Min left
                left.setMin(0.0);

                //Max left
                if (partition.getMax() > partition.getElement()->getMu() + dist)
                    left.setMax(partition.getElement()->getMu() + dist);
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(partition.getElement()->getMu() - dist);

                //Max right
                if (partition.getMax() > partition.getElement()->getCoverage() + dist)
                    right.setMax(partition.getElement()->getCoverage() + dist);
                else
                    right.setMax(partition.getMax());

            }
            else if (sqPosition == PositionRelativeToPartition::RING)
            {

                //Min left
                left.setMin(dist - partition.getElement()->getMu());

                //Max left
                if (partition.getMax() > dist + partition.getElement()->getMu())
                    left.setMax(dist + partition.getElement()->getMu());
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(0.0);

                //Max right
                if (partition.getMax() > partition.getElement()->getCoverage() + dist)
                    right.setMax(partition.getElement()->getCoverage() + dist);
                else
                    right.setMax(partition.getMax());

            }
            else
            {

                //Min left
                left.setMin(dist - partition.getElement()->getMu());

                //Max left
                if (partition.getMax() > dist + partition.getElement()->getMu())
                    left.setMax(dist + partition.getElement()->getMu());
                else
                    left.setMax(partition.getMax());

                //Min right
                right.setMin(dist - partition.getElement()->getCoverage());

                //Max right
                if (partition.getMax() > dist + partition.getElement()->getCoverage())
                    right.setMax(dist + partition.getElement()->getCoverage());
                else
                    right.setMax(partition.getMax());

            }

            return std::make_pair(left, right);

        }


    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount,iowrite,ioread";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning,iowrite,ioread";
        }

    public:
        VPTree()
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
            this->indexType = INDEX_TYPE::VPTREE_t;
            this->indexName = "VPTREE";
            this->indexFolder = "";
        }

        VPTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
               std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
               std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0,
               bool _storePivotsInLeaf = true, bool _storeDirectoryNode = false, bool _storeLeafNode = false, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->root = nullptr;
            this->pageSize = _pageSize;
            this->numPivots = _numPivots;
            this->numPerLeaf = _numPerLeaf;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->storePivotsInLeaf = _storePivotsInLeaf;
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->indexType = INDEX_TYPE::VPTREE_t;
            this->indexName = "VPTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("vp_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit VPTree(std::string _folder, std::string serializedFile = "")
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
            this->indexType = INDEX_TYPE::VPTREE_t;
            this->indexName = "VPTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~VPTree() override
        {
            deleteRecursive(std::move(root));
        }

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            size_t currentNodeID = 0;
            std::queue<std::pair<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>>> nodeQueue;
            std::pair<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>> currentNode;
            std::vector<std::pair<double, size_t>> distances;
            double dist, mu;

            if (this->dataset->getCardinality() <= numPerLeaf)
                root = std::make_unique<LeafNode<O, T>>();
            else
                root = std::make_unique<DirectoryNode<O, T>>();

            nodeQueue.push(std::make_pair(root.get(), std::move(this->dataset)));

            while (!nodeQueue.empty())
            {

                currentNode = std::move(nodeQueue.front());
                nodeQueue.pop();

                //this->pivots->clear();
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

                size_t n = distances.size();
                if (n % 2 == 0)
                    mu = (distances[n / 2 - 1].first + distances[n / 2].first) / 2.0;
                else
                    mu = distances[n / 2].first;

                if (currentNode.second->getCardinality() <= numPerLeaf)
                {

                    auto* leafNode = (LeafNode<O, T>*) currentNode.first;
                    leafNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(this->pivots->getPivot(0)));
                    leafNode->setNodeID(currentNodeID++);
                    leafNode->setMu(mu);
                    leafNode->setCoverage(distances.back().first);
                    std::filesystem::path leafIndexPath(this->indexFolder);
                    leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                    std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(this->distanceFunction->getDistanceType());
                    std::unique_ptr<pivots::Pivot<O, T>> pv = pivots::PivotFactory<O, T>::createPivot(this->pivots->getPivotType(), this->pivots);
                    std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(std::move(currentNode.second), std::move(df), std::move(pv), this->numPivots, leafIndexPath);
                    leafNode->setIndex(std::move(idx));

//                    leafNode->setDataset(std::move(currentNode.second));

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

                    std::vector<std::pair<double, size_t>> left, right;

                    for (auto& pair : distances)
                    {
                        if (pair.first <= mu)
                            left.push_back(pair);
                        else
                            right.push_back(pair);
                    }

                    auto *directoryNode = (DirectoryNode<O, T> *) currentNode.first;
                    directoryNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(this->pivots->getPivot(0)));
                    directoryNode->setNodeID(currentNodeID++);
                    directoryNode->setMu(mu);
                    directoryNode->setCoverage(distances.back().first);

                    std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>(), rightDataset = std::make_unique<dataset::Dataset<O, T>>();
                    leftDataset->setPath(currentNode.second->getPath());
                    rightDataset->setPath(currentNode.second->getPath());
                    leftDataset->setSeed(currentNode.second->getSeed());
                    rightDataset->setSeed(currentNode.second->getSeed());

                    for (auto &pair : left)
                        leftDataset->insert(currentNode.second->operator[](pair.second));

                    for (auto &pair : right)
                        rightDataset->insert(currentNode.second->operator[](pair.second));

                    if (leftDataset->getCardinality() <= numPerLeaf)
                    {
                        std::unique_ptr<vptree::Node<O, T>> leftLeafNode = std::make_unique<LeafNode<O, T>>();
                        nodeQueue.push(std::make_pair(leftLeafNode.get(), std::move(leftDataset)));
                        directoryNode->setLeft(std::move(leftLeafNode));
                    }
                    else
                    {
                        std::unique_ptr<vptree::Node<O, T>> leftDirectoryNode = std::make_unique<DirectoryNode<O, T>>();
                        nodeQueue.push(std::make_pair(leftDirectoryNode.get(), std::move(leftDataset)));
                        directoryNode->setLeft(std::move(leftDirectoryNode));
                    }

                    if (rightDataset->getCardinality() <= numPerLeaf)
                    {
                        std::unique_ptr<vptree::Node<O, T>> rightLeafNode = std::make_unique<LeafNode<O, T>>();
                        nodeQueue.push(std::make_pair(rightLeafNode.get(), std::move(rightDataset)));
                        directoryNode->setRight(std::move(rightLeafNode));
                    }
                    else
                    {
                        std::unique_ptr<vptree::Node<O, T>> rightDirectoryNode = std::make_unique<DirectoryNode<O, T>>();
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

                }

                distances.clear();

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

        void clear() override
        {
            clearRecursive(root);
        }

        std::vector<Node<O, T>*> getLeafNodes()
        {
            std::vector<Node<O, T>*> ans;
            std::queue<Node<O, T>*> nodeQueue;
            nodeQueue.push(root.get());

            while (!nodeQueue.empty())
            {

                Node<O, T>* currentNode = nodeQueue.front();
                nodeQueue.pop();

                if (currentNode->isLeafNode())
                    ans.push_back(currentNode);
                else
                {
                    nodeQueue.push(currentNode->getLeft().get());
                    nodeQueue.push(currentNode->getRight().get());
                }

            }

            return ans;
        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("VPTree::kNN not implemented yet");
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
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

                        if (!storePivotsInLeaf)
                        {
                            dist = this->distanceFunction->operator()(query, *currentLeafNode->getPivot());
                            elementQueue.push(query::ResultEntry<O>(currentLeafNode->getPivot()->getOID(), dist));
                        }

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
                        dist = this->distanceFunction->operator()(query, *currentNode->getPivot());

                        if (!storePivotsInLeaf)
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot()->getOID(), dist));

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

                        if (!storePivotsInLeaf)
                        {
                            dist = this->distanceFunction->operator()(query, *currentLeafNode->getPivot());
                            elementQueue.push(query::ResultEntry<O>(currentLeafNode->getPivot()->getOID(), dist));
                        }

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
                        dist = this->distanceFunction->operator()(query, *currentNode->getPivot());

                        if (!storePivotsInLeaf)
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot()->getOID(), dist));

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
                this->saveResultToFile(ans, query, "kNNIncremental", expt_id, {expt_id, std::to_string(k), "-1",
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

        std::unique_ptr<Node<O, T>>& getRoot()
        {
            return root;
        }

        void print(std::ostream& os) const override
        {

            std::stack<std::pair<Node<O, T>*, size_t>> nodeStack;
            nodeStack.push(std::make_pair(root.get(), 0));

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "VPTree" << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Number of objects per leaf: " << numPerLeaf << std::endl;
            os << "Store pivots in leaf: " << (storePivotsInLeaf ? "true" : "false") << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;

            while (!nodeStack.empty())
            {

                auto currentNode = nodeStack.top();
                nodeStack.pop();

                os << "**********************************************************************************************************************************************************************************\n\n";
                os << *currentNode.first << std::endl;

                if (currentNode.first->getLeft() != nullptr)
                    nodeStack.push(std::make_pair(currentNode.first->getLeft().get(), currentNode.second + 1));

                if (currentNode.first->getRight() != nullptr)
                    nodeStack.push(std::make_pair(currentNode.first->getRight().get(), currentNode.second + 1));

            }

        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<VPTree<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

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
            ans += sizeof(size_t) * 5;

            ans += sizeof(size_t) + serializedTree.size();

            return ans;
        }

    };

}

#endif //GERVLIB_VPTREE_H
