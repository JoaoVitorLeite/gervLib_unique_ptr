//
// Created by joaoleite on 8/29/23.
//

#ifndef GERVLIB_KDTREE_H
#define GERVLIB_KDTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "NAryTree.h"

namespace gervLib::index::kdtree
{

    const double MIN_KDTREE = 0.0;

    template <typename O, typename T>
    class Node : public serialize::Serialize
    {

    protected:
        std::unique_ptr<Node<O, T>> left, right;
        size_t nodeID{};
        std::unique_ptr<std::vector<std::pair<double, double>>> bounds;
        MEMORY_STATUS memoryStatus = MEMORY_STATUS::NONE;

    public:
        Node()
        {
            left = nullptr;
            right = nullptr;
            nodeID = 0;
            bounds = nullptr;
        }

        Node(std::unique_ptr<Node<O, T>> _left, std::unique_ptr<Node<O, T>> _right)
        {
            this->left = std::move(_left);
            this->right = std::move(_right);
            this->nodeID = 0;
            this->bounds = nullptr;
        }

        virtual ~Node()
        {
            if (bounds != nullptr) {
                bounds->clear();
                bounds.reset();
            }
        }

        void setBoundsSize(size_t sz)
        {
            if (bounds != nullptr)
            {
                bounds->clear();
                bounds.reset();
            }

            bounds = std::make_unique<std::vector<std::pair<double, double>>>(sz, std::make_pair(MIN_KDTREE, std::numeric_limits<double>::max()));
        }


        bool equalBounds(std::unique_ptr<std::vector<std::pair<double, double>>>& other)
        {

            if ((bounds == nullptr && other != nullptr) || (bounds != nullptr && other == nullptr))
                return false;

            if (bounds == nullptr && other == nullptr)
                return true;
            else if (bounds->size() != other->size())
                return false;
            else {

                for (size_t i = 0; i < bounds->size(); i++) {
                    if (bounds->at(i).first != other->at(i).first || bounds->at(i).second != other->at(i).second)
                        return false;
                }
            }

            return true;

        }

        std::unique_ptr<std::vector<std::pair<double, double>>>& getBoundary()
        {
            return bounds;
        }

        size_t getBoundsSize()
        {
            if (bounds == nullptr)
                return 0;
            else
                return bounds->size();
        }

        std::unique_ptr<Node<O, T>>& getLeft()
        {
            return left;
        }

        std::unique_ptr<Node<O, T>>& getRight()
        {
            return right;
        }

        void setLeft(std::unique_ptr<Node<O, T>> _left)
        {
            this->left.reset();
            left = std::move(_left);
        }

        void setRight(std::unique_ptr<Node<O, T>> _right)
        {
            this->right.reset();
            right = std::move(_right);
        }

        void setBoundary(std::unique_ptr<std::vector<std::pair<double, double>>> _bounds)
        {
            if (bounds != nullptr) {
                bounds->clear();
                bounds.reset();
            }

            bounds = std::move(_bounds);
        }

        void setBoundary(std::vector<std::pair<double, double>>& _bounds)
        {
            if (bounds != nullptr) {
                bounds->clear();
                bounds.reset();
            }

            bounds = std::make_unique<std::vector<std::pair<double, double>>>(_bounds);
        }

        void setNodeID(size_t id)
        {
            nodeID = id;
        }

        size_t getNodeID()
        {
            return nodeID;
        }

        void setMemoryStatus(MEMORY_STATUS status)
        {
            memoryStatus = status;
        }

        MEMORY_STATUS getMemoryStatus()
        {
            return memoryStatus;
        }

        void setMinBound(size_t index, double value)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::setMinBound: index out of range");
            bounds->at(index).first = value;
        }

        void setMaxBound(size_t index, double value)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::setMaxBound: index out of range");
            bounds->at(index).second = value;
        }

        double getMinBound(size_t index)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::getMinBound: index out of range");
            return bounds->at(index).first;
        }

        double getMaxBound(size_t index)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::getMaxBound: index out of range");
            return bounds->at(index).second;
        }

        void setBound(size_t index, double min, double max)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::setBound: index out of range");
            bounds->at(index).first = min;
            bounds->at(index).second = max;
        }

        std::pair<double, double> getBound(size_t index)
        {
            utils::check_range(0, bounds->size()-1, index, "Node::getBound: index out of range");
            return bounds->at(index);
        }

        virtual bool isLeafNode()
        {
            return false;
        }

        virtual bool isEqual(std::unique_ptr<Node<O, T>>& other)
        {
            if (memoryStatus == MEMORY_STATUS::IN_DISK || other->memoryStatus == MEMORY_STATUS::IN_DISK)
                return nodeID == other->nodeID;
            else
            {
                return equalBounds(other->bounds);
            }
        }

        virtual void clear()
        {
            if (bounds != nullptr)
            {
                bounds->clear();
                bounds.reset();
            }
        }

        virtual void print(std::ostream& os) const
        {

            os << "Node type: Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[memoryStatus] << std::endl;
            os << "Bounds: " << std::endl;

            if (bounds != nullptr)
            {
                for (size_t i = 0; i < bounds->size(); i++)
                    os << "Dimension " << i << ": [" << bounds->at(i).first << ", " << bounds->at(i).second << "]" << std::endl;
            }
            else
                os << "Bounds not set" << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux = memoryStatusMap[memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;

            sz = (bounds == nullptr ? 0 : bounds->size());

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            if (bounds != nullptr)
            {
                for (size_t i = 0; i < bounds->size(); i++)
                {
                    memcpy(data.get() + offset, &bounds->at(i).first, sizeof(double));
                    offset += sizeof(double);
                    memcpy(data.get() + offset, &bounds->at(i).second, sizeof(double));
                    offset += sizeof(double);
                }
            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            memoryStatus = memoryStatusMapReverse[aux];

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                if (bounds != nullptr)
                {
                    bounds->clear();
                    bounds.reset();
                }

                bounds = std::make_unique<std::vector<std::pair<double, double>>>(sz, std::make_pair<double, double>(0.0, 0.0));

                for (size_t i = 0; i < sz; i++)
                {
                    memcpy(&bounds->at(i).first, _data.get() + offset, sizeof(double));
                    offset += sizeof(double);
                    memcpy(&bounds->at(i).second, _data.get() + offset, sizeof(double));
                    offset += sizeof(double);
                }
            }
            else
                bounds = nullptr;

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;
            ans += sizeof(size_t); //nodeID
            ans += sizeof(size_t) + memoryStatusMap[memoryStatus].size();
            ans += sizeof(size_t) + (bounds == nullptr ? 0 : sizeof(double) * 2 * bounds->size());

            return ans;

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

        ~DirectoryNode() override = default;

        bool isLeafNode() override
        {
            return false;
        }

        void print(std::ostream& os) const override
        {
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;
            os << "Bounds: " << std::endl;

            if (this->bounds != nullptr)
            {
                for (size_t i = 0; i < this->bounds->size(); i++)
                    os << "Dimension " << i << ": [" << this->bounds->at(i).first << ", " << this->bounds->at(i).second << "]" << std::endl;
            }
            else
                os << "Bounds not set" << std::endl;
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

            if (!this->equalBounds(node->getBoundary()))
                return false;

            return true;

        }

        void print(std::ostream& os) const override
        {

            os << "Node type: Leaf Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;
            os << "Bounds: " << std::endl;

            if (this->bounds != nullptr)
            {
                for (size_t i = 0; i < this->bounds->size(); i++)
                    os << "Dimension " << i << ": [" << this->bounds->at(i).first << ", " << this->bounds->at(i).second << "]" << std::endl;
            }
            else
                os << "Bounds not set" << std::endl;

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
    class KdTree : public Index<O, T>
    {

    private:
        std::unique_ptr<Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
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

        double medianByDimension(std::unique_ptr<dataset::Dataset<O, T>>& dataset, size_t dim)
        {
            std::vector<double> vec(dataset->getCardinality());

            for (size_t i = 0; i < dataset->getCardinality(); i++)
                vec[i] = dataset->getElement(i).operator[](dim);

            std::sort(vec.begin(), vec.end());

            if(vec.size() % 2 == 0)
            {

                return (vec[vec.size()/2 - 1] + vec[vec.size()/2])/2.0;

            }
            else
            {

                return vec[vec.size()/2];

            }

        }

        std::pair<std::unique_ptr<std::vector<std::pair<double, double>>>, std::unique_ptr<std::vector<std::pair<double, double>>>> splitBounds(std::unique_ptr<std::vector<std::pair<double, double>>>& bound,
                                                                                                                                                double median, size_t pos)
        {
            std::unique_ptr<std::vector<std::pair<double, double>>> cp1 = std::make_unique<std::vector<std::pair<double, double>>>(*bound);
            std::unique_ptr<std::vector<std::pair<double, double>>> cp2 = std::make_unique<std::vector<std::pair<double, double>>>(*bound);
            std::pair<double, double> aux_cp1 = std::make_pair(bound->at(pos).first, std::nextafter(median, -std::numeric_limits<double>::infinity())),
                                      aux_cp2 = std::make_pair(median, bound->at(pos).second);
            cp1->at(pos) = aux_cp1;
            cp2->at(pos) = aux_cp2;
            return std::make_pair(std::move(cp1), std::move(cp2));

        }
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

        bool isInterval(double infBound, double supBound, double test)
        {

            return ((test >= infBound) && (test <= supBound));

        }

        double minDist(dataset::BasicArrayObject<O, T>& query, std::unique_ptr<std::vector<std::pair<double, double>>>& bounds)
        {
            double limInfCase3 = -1.0;
            double limInfCase2 = -1.0;
            double answer = -1.0;
            bool within = true;

            for(size_t x = 0; x < bounds->size(); x++)
            {

                if(!isInterval(std::abs(bounds->at(x).first), std::abs(bounds->at(x).second), query.operator[](x)))
                {

                    within = false;

                    limInfCase3 = std::max(limInfCase3,

                                           std::min(

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).first)),

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).second))

                                           )

                    );
                }
                else
                {

                    limInfCase2 = std::min(limInfCase2,

                                           std::min(

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).first)),

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).second))

                                           )

                    );

                }

            }

            if(within)
            {

                answer = 0.0;

            }
            else
            {

                if(limInfCase2 != -1.0)
                {

                    answer = limInfCase2;

                }
                else
                {

                    answer = limInfCase3;

                }

            }

            return answer;
        }

        double maxDist(dataset::BasicArrayObject<O, T>& query, std::unique_ptr<std::vector<std::pair<double, double>>>& bounds)
        {
            double answer = -1.0;

            answer = std::numeric_limits<double>::max();

            for(size_t x = 0; x < bounds->size(); x++)
            {

                if (std::numeric_limits<double>::max() -  std::abs(bounds->at(x).first) >= query.operator[](x))
                {

                    answer = std::min(answer, query.operator[](x) + std::abs(bounds->at(x).first));

                }

                if (std::numeric_limits<double>::max() -  std::abs(bounds->at(x).second) >= query.operator[](x))
                {

                    answer = std::min(answer, query.operator[](x) + std::abs(bounds->at(x).second));

                }

            }

            return answer;
        }

    public:
        KdTree()
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
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::KDTREE;
            this->indexName = "KDTREE";
            this->indexFolder = "";
        }

        KdTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
               std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
               std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0,
               bool _storeDirectoryNode = false, bool _storeLeafNode = false, bool _useLAESA = false, std::string folder="")
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
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::KDTREE;
            this->indexName = "KDTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("kd_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit KdTree(std::string _folder, std::string serializedFile = "")
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
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::KDTREE;
            this->indexName = "KDTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~KdTree() override
        {
            deleteRecursive(std::move(root));
        }

        std::unique_ptr<Node<O, T>>& getRoot()
        {
            return root;
        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<KdTree<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

        }

        void print(std::ostream& os) const override {

            std::stack<std::pair<Node<O, T> *, size_t>> nodeStack;
            nodeStack.push(std::make_pair(root.get(), 0));

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "KdTree" << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Number of objects per leaf: " << numPerLeaf << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;
            os << "Use LAESA: " << (useLAESA ? "true" : "false") << std::endl;

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

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            size_t currentNodeID = 0;
            std::queue<std::tuple<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>, size_t>> nodeQueue;
            std::tuple<Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>, size_t> currentTuple;
            Node<O, T>* currentNode;
            std::unique_ptr<dataset::Dataset<O, T>> currentDataset;
            size_t dPartition;
            double median;
            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);

            if (this->dataset->getCardinality() <= numPerLeaf)
                root = std::make_unique<LeafNode<O, T>>();
            else
                root = std::make_unique<DirectoryNode<O, T>>();

            root->setBoundsSize(this->numPivots);
            root->setNodeID(currentNodeID++);

            nodeQueue.push(std::make_tuple(root.get(), std::move(this->dataset), 0));

            while (!nodeQueue.empty()) {

                currentTuple = std::move(nodeQueue.front());
                nodeQueue.pop();

                currentNode = std::get<0>(currentTuple);
                currentDataset = std::move(std::get<1>(currentTuple));
                dPartition = std::get<2>(currentTuple);

                if (currentNode->isLeafNode())
                {

                    auto *leafNode = dynamic_cast<LeafNode<O, T>*>(currentNode);

                    if (useLAESA) {
                        std::filesystem::path leafIndexPath(this->indexFolder);
                        leafIndexPath /= "laesa_leafnode_" + std::to_string(currentNode->getNodeID());
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                this->distanceFunction->getDistanceType());
                        std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                std::move(currentDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                this->numPivots,
                                leafIndexPath);
                        this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + idx->getDistanceCount());
                        leafNode->setIndex(std::move(idx));
                    }
                    else
                    {
                        leafNode->setDataset(std::move(currentDataset));
                    }
                }
                else
                {
                    median = medianByDimension(currentDataset, dPartition);

                    std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>(), rightDataset = std::make_unique<dataset::Dataset<O, T>>();
                    leftDataset->setSeed(currentDataset->getSeed());
                    rightDataset->setSeed(currentDataset->getSeed());
                    leftDataset->setPath(currentDataset->getPath());
                    rightDataset->setPath(currentDataset->getPath());
                    leftDataset->setDimensionality(currentDataset->getDimensionality());
                    rightDataset->setDimensionality(currentDataset->getDimensionality());

                    for (size_t i = 0; i < currentDataset->getCardinality(); i++) {

                        if (currentDataset->getElement(i).operator[](dPartition) <= median)
                            leftDataset->insert(currentDataset->getElement(i));
                        else
                            rightDataset->insert(currentDataset->getElement(i));

                    }

                    std::pair<std::unique_ptr<std::vector<std::pair<double, double>>>, std::unique_ptr<std::vector<std::pair<double, double>>>> splitBoundaries = splitBounds(currentNode->getBoundary(), median, dPartition);

                    if (leftDataset->getCardinality() <= numPerLeaf)
                    {

                        std::unique_ptr<LeafNode<O, T>> leafNode = std::make_unique<LeafNode<O, T>>();
                        leafNode->setBoundary(std::move(splitBoundaries.first));
                        leafNode->setNodeID(currentNodeID++);

                        if (useLAESA) {
                            std::filesystem::path leafIndexPath(this->indexFolder);
                            leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                    this->distanceFunction->getDistanceType());
                            std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                    std::move(leftDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                    this->numPivots,
                                    leafIndexPath);
                            this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + idx->getDistanceCount());
                            leafNode->setIndex(std::move(idx));
                        }
                        else
                        {
                            leafNode->setDataset(std::move(leftDataset));
                        }

                        if (storeLeafNode)
                        {
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> leafData = leafNode->serialize();
                            this->pageManager->save(leafNode->getNodeID(), std::move(leafData), leafNode->getSerializedSize());
                            leafNode->clear();
                        }
                        else
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setLeft(std::move(leafNode));

                    }
                    else
                    {

                        std::unique_ptr<DirectoryNode<O, T>> directoryNode = std::make_unique<DirectoryNode<O, T>>();
                        directoryNode->setBoundary(std::move(splitBoundaries.first));
                        directoryNode->setNodeID(currentNodeID++);

                        if (this->storeDirectoryNode) {
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                            this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                            directoryNode->clear();
                        }
                        else
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setLeft(std::move(directoryNode));
                        nodeQueue.push(std::make_tuple(currentNode->getLeft().get(), std::move(leftDataset), (dPartition + 1) % currentDataset->getDimensionality()));

                    }

                    if (rightDataset->getCardinality() <= numPerLeaf)
                    {

                        std::unique_ptr<LeafNode<O, T>> leafNode = std::make_unique<LeafNode<O, T>>();
                        leafNode->setBoundary(std::move(splitBoundaries.second));
                        leafNode->setNodeID(currentNodeID++);

                        if (useLAESA) {
                            std::filesystem::path leafIndexPath(this->indexFolder);
                            leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                    this->distanceFunction->getDistanceType());
                            std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                    std::move(rightDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                    this->numPivots, leafIndexPath);
                            this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + idx->getDistanceCount());
                            leafNode->setIndex(std::move(idx));
                        }
                        else
                        {
                            leafNode->setDataset(std::move(rightDataset));
                        }

                        if (storeLeafNode)
                        {
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> leafData = leafNode->serialize();
                            this->pageManager->save(leafNode->getNodeID(), std::move(leafData), leafNode->getSerializedSize());
                            leafNode->clear();
                        }
                        else
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setRight(std::move(leafNode));

                    }
                    else
                    {

                        std::unique_ptr<DirectoryNode<O, T>> directoryNode = std::make_unique<DirectoryNode<O, T>>();
                        directoryNode->setBoundary(std::move(splitBoundaries.second));
                        directoryNode->setNodeID(currentNodeID++);

                        if (this->storeDirectoryNode) {
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                            this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                            directoryNode->clear();
                        }
                        else
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setRight(std::move(directoryNode));
                        nodeQueue.push(std::make_tuple(currentNode->getRight().get(), std::move(rightDataset), (dPartition + 1) % currentDataset->getDimensionality()));

                    }

                }

                if (currentNode->isLeafNode())
                {
                    if (storeLeafNode)
                    {
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> leafData = currentNode->serialize();
                        this->pageManager->save(currentNode->getNodeID(), std::move(leafData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);
                }
                else
                {
                    if (this->storeDirectoryNode) {
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> directoryData = currentNode->serialize();
                        this->pageManager->save(currentNode->getNodeID(), std::move(directoryData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);
                }

                currentDataset->clear();
                currentDataset.reset();

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

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults, bool saveStatistics) override
        {

            utils::Timer timer{};
            timer.start();
            std::string expt_id = utils::generateExperimentID();
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

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty())) {

                if (elementQueue.empty())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->isLeafNode())
                    {

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                            currentNode->deserialize(std::move(nodeData));
                        }

                        currentLeafNode = (LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k, true);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();
                            this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + currentLeafNode->getIndex()->getDistanceFunction()->getDistanceCount());

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->clear();
                        }

                    }
                    else
                    {

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getLeft()->getNodeID());
                            currentNode->getLeft()->deserialize(std::move(nodeData));
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getRight()->getNodeID());
                            currentNode->getRight()->deserialize(std::move(nodeData));
                        }

                        nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getLeft().get(), minDist(query, currentNode->getLeft()->getBoundary()), maxDist(query, currentNode->getLeft()->getBoundary())));
                        nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getRight().get(), minDist(query, currentNode->getRight()->getBoundary()), maxDist(query, currentNode->getRight()->getBoundary())));

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getLeft()->clear();
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getRight()->clear();
                        }

                    }

                }
                else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->isLeafNode())
                    {

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                            currentNode->deserialize(std::move(nodeData));
                        }

                        currentLeafNode = (LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k, true);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();
                            this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + currentLeafNode->getIndex()->getDistanceFunction()->getDistanceCount());

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->clear();
                        }

                    }
                    else
                    {

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getLeft()->getNodeID());
                            currentNode->getLeft()->deserialize(std::move(nodeData));
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getRight()->getNodeID());
                            currentNode->getRight()->deserialize(std::move(nodeData));
                        }

                        nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getLeft().get(), minDist(query, currentNode->getLeft()->getBoundary()), maxDist(query, currentNode->getLeft()->getBoundary())));
                        nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getRight().get(), minDist(query, currentNode->getRight()->getBoundary()), maxDist(query, currentNode->getRight()->getBoundary())));

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getLeft()->clear();
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getRight()->clear();
                        }

                    }

                }
                else
                {

                    result.push(elementQueue.top());
                    elementQueue.pop();

                    if (result.size() % 5 == 0 && result.size() != 100)
                    {
                        timer.stop();

                        if (saveStatistics)
                        {
                            this->saveStatistics({expt_id, std::to_string(result.size()), "-1",
                                                  std::to_string(timer.getElapsedTime()),
                                                  std::to_string(timer.getElapsedTimeSystem()),
                                                  std::to_string(timer.getElapsedTimeUser()),
                                                  std::to_string(this->distanceFunction->getDistanceCount()),
                                                  std::to_string(this->prunning),
                                                  std::to_string(configure::IOWrite - ioW),
                                                  std::to_string(configure::IORead - ioR)});
                        }

                    }

                }

            }

            std::vector<query::ResultEntry<O>> ans = result.getResults();
            std::reverse(ans.begin(), ans.end());

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
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
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

            sz = (storeDirectoryNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storeLeafNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (useLAESA ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = sizeof(size_t) + serializedTree.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, serializedTree.c_str(), serializedTree.size());
            offset += serializedTree.size();

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
            ans += sizeof(size_t) * 5;

            ans += sizeof(size_t) + serializedTree.size();

            return ans;
        }

    };

}

#endif //GERVLIB_KDTREE_H
