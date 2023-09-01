//
// Created by joaoleite on 9/1/23.
//

#ifndef GERVLIB_PMTREE_H
#define GERVLIB_PMTREE_H

#include "Index.h"
#include "IndexFactory.h"

namespace gervLib::index::pmtree
{

    template <typename O, typename T>
    class Node : public serialize::Serialize
    {

    protected:
        Node<O, T>* parent;
        std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot;
        double dist_to_parent, coverage;
        std::unique_ptr<std::vector<std::pair<double, double>>> rings;
        size_t nodeID;
        std::vector<std::unique_ptr<Node<O, T>>> childrens;
        MEMORY_STATUS memoryStatus;

    public:
        Node()
        {
            parent = nullptr;
            pivot = nullptr;
            dist_to_parent = 0.0;
            coverage = 0.0;
            nodeID = 0;
            rings = nullptr;
            memoryStatus = MEMORY_STATUS::IN_MEMORY;
        }

        explicit Node(size_t num_pivots)
        {
            parent = nullptr;
            pivot = nullptr;
            dist_to_parent = 0.0;
            coverage = 0.0;
            nodeID = 0;
            rings = std::make_unique<std::vector<std::pair<double, double>>>(num_pivots, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
            memoryStatus = MEMORY_STATUS::IN_MEMORY;
        }

        void setNumberOfPivots(size_t num_pivots)
        {

            if (rings != nullptr)
            {
                rings->clear();
                rings.reset();
            }

            rings = std::make_unique<std::vector<std::pair<double, double>>>(num_pivots, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

        }

        void setParent(Node<O, T>* _parent)
        {
            this->parent = _parent;
        }

        Node<O, T>* getParent()
        {
            return parent;
        }

        virtual ~Node()
        {
            if (this->pivot != nullptr) {
                this->pivot->clear();
                this->pivot.reset();
            }

            if (this->rings != nullptr) {
                this->rings->clear();
                this->rings.reset();
            }
        }

        virtual void clear()
        {
            if (this->pivot != nullptr) {
                this->pivot->clear();
                this->pivot.reset();
            }

            if (this->rings != nullptr) {
                this->rings->clear();
                this->rings.reset();
            }
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

        std::unique_ptr<dataset::BasicArrayObject<O, T>>& getPivot()
        {
            return pivot;
        }

        void setDistanceToParent(double _dist_to_parent)
        {
            this->dist_to_parent = _dist_to_parent;
        }

        double getDistanceToParent()
        {
            return dist_to_parent;
        }

        void setCoverage(double _coverage)
        {
            this->coverage = _coverage;
        }

        double getCoverage()
        {
            return coverage;
        }

        void setNodeID(size_t _nodeID)
        {
            this->nodeID = _nodeID;
        }

        size_t getNodeID()
        {
            return nodeID;
        }

        void setRing(size_t index, double min, double max)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            rings->at(index) = std::make_pair(min, max);
        }

        void setRingMin(size_t index, double min)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            rings->at(index).first = min;
        }

        void setRingMax(size_t index, double max)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            rings->at(index).second = max;
        }

        double getRingMin(size_t index)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            return rings->at(index).first;
        }

        double getRingMax(size_t index)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            return rings->at(index).second;
        }

        std::pair<double, double> getRing(size_t index)
        {
            utils::check_range(0, rings->size()-1, index, "Index out of bounds");
            return rings->at(index);
        }

        void setChild(size_t index, std::unique_ptr<Node<O, T>> child)
        {
            utils::check_range(0, childrens.size()-1, index, "Index out of bounds");
            childrens[index] = std::move(child);
        }

        void addChild(std::unique_ptr<Node<O, T>> child)
        {
            child->setParent(this);
            childrens.push_back(std::move(child));
        }

        std::unique_ptr<Node<O, T>>& getChild(size_t index)
        {
            utils::check_range(0, childrens.size()-1, index, "Index out of bounds");
            return childrens[index];
        }

        std::vector<std::unique_ptr<Node<O, T>>>& getChildrens()
        {
            return childrens;
        }

        size_t getChildrensSize()
        {
            if (childrens.empty())
                return 0;
            else
                return childrens.size();
        }

        void setMemoryStatus(MEMORY_STATUS _memoryStatus)
        {
            this->memoryStatus = _memoryStatus;
        }

        MEMORY_STATUS getMemoryStatus()
        {
            return memoryStatus;
        }

        virtual bool isEqual(std::unique_ptr<Node<O, T>> &other)
        {
            if (memoryStatus == index::MEMORY_STATUS::IN_DISK)
            {
                return nodeID == other->nodeID;
            }
            else
            {

                if ((pivot == nullptr && other->pivot != nullptr) || (pivot != nullptr && other->pivot == nullptr))
                    return false;
                else if (pivot != nullptr && other->pivot != nullptr)
                    if (!pivot->isEqual(*other->pivot))
                        return false;

                if ((rings == nullptr && other->rings != nullptr) || (rings != nullptr && other->rings == nullptr))
                    return false;
                else if (rings != nullptr && other->rings != nullptr)
                    if (*rings != *other->rings)
                        return false;

                if (dist_to_parent != other->dist_to_parent)
                    return false;

                if (coverage != other->coverage)
                    return false;

                return true;

            }

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
            os << "Node type: Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;
            os << "Distance to parent: " << dist_to_parent << std::endl;
            os << "Coverage: " << coverage << std::endl;

            if (pivot != nullptr)
                os << "Pivot: " << *pivot << std::endl;
            else
                os << "Pivot: Null" << std::endl;

            if (rings != nullptr)
            {
                os << "Rings: ";
                for (size_t i = 0; i < rings->size(); i++)
                    os << "(" << rings->at(i).first << ", " << rings->at(i).second << ")\n";
                os << std::endl;
            }
            else
                os << "Rings: Null" << std::endl;

            os << "Number of Childrens: " << childrens.size() << std::endl;
            os << "Memory status: " << index::memoryStatusMap[memoryStatus] << std::endl;
        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            if (pivot != nullptr)
            {
                sz = pivot->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> pivot_data = pivot->serialize();
                memcpy(data.get() + offset, pivot_data.get(), pivot->getSerializedSize());
                offset += pivot->getSerializedSize();
                pivot_data.reset();
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            memcpy(data.get() + offset, &dist_to_parent, sizeof(double));
            offset += sizeof(double);

            memcpy(data.get() + offset, &coverage, sizeof(double));
            offset += sizeof(double);

            std::string aux = index::memoryStatusMap[memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), aux.size());
            offset += aux.size();
            aux.clear();

            if (rings != nullptr)
            {
                sz = rings->size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                for (size_t i = 0; i < rings->size(); i++)
                {
                    memcpy(data.get() + offset, &rings->at(i).first, sizeof(double));
                    offset += sizeof(double);

                    memcpy(data.get() + offset, &rings->at(i).second, sizeof(double));
                    offset += sizeof(double);
                }
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

            memcpy(&nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                std::unique_ptr<u_char[]> pivot_data = std::make_unique<u_char[]>(sz);
                memcpy(pivot_data.get(), _data.get() + offset, sz);
                offset += sz;

                if (pivot != nullptr)
                {
                    pivot->clear();
                    pivot.reset();
                }

                pivot = std::make_unique<dataset::BasicArrayObject<O, T>>();
                pivot->deserialize(std::move(pivot_data));
            }
            else
                pivot = nullptr;

            memcpy(&dist_to_parent, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&coverage, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            memoryStatus = index::memoryStatusMapReverse[aux];

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {

                if (rings != nullptr)
                {
                    rings->clear();
                    rings.reset();
                }

                rings = std::make_unique<std::vector<std::pair<double, double>>>(sz);

                for (size_t i = 0; i < sz; i++)
                {
                    memcpy(&rings->at(i).first, _data.get() + offset, sizeof(double));
                    offset += sizeof(double);

                    memcpy(&rings->at(i).second, _data.get() + offset, sizeof(double));
                    offset += sizeof(double);
                }
            }
            else
                rings = nullptr;

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = sizeof(size_t); //nodeID
            ans += sizeof(size_t) + (pivot != nullptr ? pivot->getSerializedSize() : 0); //pivot
            ans += sizeof(double) * 2; // dist_to_parent, coverage
            ans += sizeof(size_t) + memoryStatusMap[memoryStatus].size(); //memoryStatus
            ans += sizeof(size_t) + (rings != nullptr ? rings->size() * sizeof(double) * 2 : 0); //rings

            return ans;
        }

    };

    template <typename O, typename T>
    class DirectoryNode : public Node<O, T>
    {

    public:
        explicit DirectoryNode(size_t num_pivots) : Node<O, T>(num_pivots) { }

        DirectoryNode() : Node<O, T>() { }

        ~DirectoryNode() override = default;

        bool isLeafNode() override
        {
            return false;
        }

        void print(std::ostream& os) const
        {
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Distance to parent: " << this->dist_to_parent << std::endl;
            os << "Coverage: " << this->coverage << std::endl;

            if (this->pivot != nullptr)
                os << "Pivot: " << *this->pivot << std::endl;
            else
                os << "Pivot: Null" << std::endl;

            if (this->rings != nullptr)
            {
                os << "Rings: ";
                for (size_t i = 0; i < this->rings->size(); i++)
                    os << "(" << this->rings->at(i).first << ", " << this->rings->at(i).second << ")\n";
                os << std::endl;
            }
            else
                os << "Rings: Null" << std::endl;

            os << "Number of Childrens: " << this->childrens.size() << std::endl;
            os << "Memory status: " << index::memoryStatusMap[this->memoryStatus] << std::endl;
        }

    };

    template <typename O, typename T>
    class LeafNode : public Node<O, T>
    {

    private:
        std::unique_ptr<Index<O, T>> index;
        std::unique_ptr<dataset::Dataset<O, T>> dataset;

    public:
        LeafNode() : Node<O, T>()
        {
            index = nullptr;
            dataset = nullptr;
        }

        explicit LeafNode(size_t num_pivots) : Node<O, T>(num_pivots)
        {
            index = nullptr;
            dataset = nullptr;
        }

        LeafNode(size_t num_pivots, std::unique_ptr<Index<O, T>> _index) : Node<O, T>(num_pivots)
        {
            index = std::move(_index);
            dataset = nullptr;
        }

        LeafNode(size_t num_pivots, std::unique_ptr<dataset::Dataset<O, T>> _dataset) : Node<O, T>(num_pivots)
        {
            index = nullptr;
            dataset = std::move(_dataset);
        }

        ~LeafNode() override
        {
            Node<O, T>::clear();
            if (this->dataset != nullptr) {
                this->dataset->clear();
                this->dataset.reset();
            }
            if (this->index != nullptr) {
                this->index->clear();
                this->index.reset();
            }
        }

        bool isLeafNode() override
        {
            return true;
        }

        size_t getLeafSize()
        {
            if (dataset == nullptr)
                return 0;
            else
                return dataset->getCardinality();
        }

        void clear() override
        {
            Node<O, T>::clear();
            if (this->dataset != nullptr) {
                this->dataset->clear();
                this->dataset.reset();
            }
            if (this->index != nullptr) {
                this->index->clear();
                this->index.reset();
            }
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

        void insert(dataset::BasicArrayObject<O, T> object)
        {
            if (dataset == nullptr)
                dataset = std::make_unique<dataset::Dataset<O, T>>();
            this->dataset->insert(object);
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

            if ((this->rings == nullptr && node->rings != nullptr) || (this->rings != nullptr && node->rings == nullptr))
                return false;
            else if (this->rings != nullptr && node->rings != nullptr)
                if (*this->rings != *node->rings)
                    return false;

            if (this->dist_to_parent != node->dist_to_parent)
                return false;

            if (this->coverage != node->coverage)
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
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
            os << "Distance to parent: " << this->dist_to_parent << std::endl;
            os << "Coverage: " << this->coverage << std::endl;

            if (this->pivot != nullptr)
                os << "Pivot: " << *this->pivot << std::endl;
            else
                os << "Pivot: Null" << std::endl;

            if (this->rings != nullptr)
            {
                os << "Rings: ";
                for (size_t i = 0; i < this->rings->size(); i++)
                    os << "(" << this->rings->at(i).first << ", " << this->rings->at(i).second << ")\n";
                os << std::endl;
            }
            else
                os << "Rings: Null" << std::endl;

            if (this->dataset != nullptr)
                os << "Dataset: " << *this->dataset << std::endl;
            else
                os << "Dataset: Null" << std::endl;

            if (this->index != nullptr)
                os << "Index: " << *this->index << std::endl;
            else
                os << "Index: Null" << std::endl;

            os << "Number of Childrens: " << this->childrens.size() << std::endl;
            os << "Memory status: " << index::memoryStatusMap[this->memoryStatus] << std::endl;
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
    class PMTree : public Index<O, T>
    {

    private:
        std::unique_ptr<Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};

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
        PMTree()
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
            this->indexType = INDEX_TYPE::PMTREE;
            this->indexName = "PMTREE";
            this->indexFolder = "";
        }

        PMTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset, std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
               std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0, bool _storeDirectoryNode = false,
               bool _storeLeafNode = false, bool _useLAESA = true, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->root = nullptr;
            this->pageSize = _pageSize;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numPerLeaf = _numPerLeaf;
            this->numPivots = _numPivots;
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::PMTREE;
            this->indexName = "PMTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("vp_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit PMTree(std::string _folder, std::string serializedFile = "")
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
            this->indexType = INDEX_TYPE::PMTREE;
            this->indexName = "PMTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);
        }

        ~PMTree() override = default;

        //TODO implement delete, clear, print, isEqual, buildIndex, kNN, kNNIncremental, serialize, deserialize, getSerializedSize

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;

            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);

            if (root != nullptr)
            {
                root->clear();
                root.reset();
            }

            root = std::make_unique<LeafNode<O, T>>(this->numPivots);

            for (size_t i = 0; i < this->dataset->getCardinality(); i++)
            {
                insert(this->dataset->getElement(i), root.get());
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

        void insert(dataset::BasicArrayObject<O, T>& element, Node<O, T>* node)
        {

            if (node->isLeafNode())
            {

                auto* leafNode = dynamic_cast<LeafNode<O, T>*>(node);

                if (leafNode->getLeafSize() >= this->numPerLeaf)
                    split(element, node);
                else
                {

                    if (leafNode->getPivot() == nullptr)
                        leafNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(element));

                    leafNode->insert(element);

                    if (leafNode->getParent() == nullptr)
                        leafNode->setDistanceToParent(std::numeric_limits<double>::max());
                    else
                        leafNode->setDistanceToParent(this->distanceFunction->operator()(*leafNode->getPivot(), *leafNode->getParent()->getPivot()));

                    leafNode->setCoverage(std::max(this->distanceFunction->operator()(*leafNode->getPivot(), element), leafNode->getCoverage()));
                    merge_rings_and_covarage(leafNode->getParent());

                }

            }
            else
            {

                std::vector<double> withIncrease, withoutIncrease;
                double dist;

                for (size_t i = 0; i < node->getChildrensSize(); i++)
                {

                    dist = this->distanceFunction->operator()(element, *node->getChild(i)->getPivot());

                    if (dist <= node->getChild(i)->getCoverage())
                        withoutIncrease.push_back(dist);
                    else
                        withIncrease.push_back(dist - node->getChild(i)->getCoverage());

                }

                if (withoutIncrease.empty())
                {
                    auto min_pos = std::min_element(withIncrease.begin(), withIncrease.end());
                    insert(element, node->getChild(std::distance(withIncrease.begin(), min_pos)).get());
                }
                else
                {
                    auto min_pos = std::min_element(withoutIncrease.begin(), withoutIncrease.end());
                    insert(element, node->getChild(std::distance(withoutIncrease.begin(), min_pos)).get());
                }

                withoutIncrease.clear();
                withIncrease.clear();

            }

        }

        void update_rings(Node<O, T>* node)
        {
            if (!node->isLeafNode()) {

                for (size_t i = 0; i < node->getChildrensSize(); i++) {

                    for (size_t p = 0; p < this->numPivots; p++) {

                        node->setRingMin(p, std::min(node->getRingMin(p), node->getChild(i)->getRingMin(p)));
                        node->setRingMax(p, std::max(node->getRingMax(p), node->getChild(i)->getRingMax(p)));

                    }

                }

            }
            else
            {
                auto* leaf = dynamic_cast<LeafNode<O, T>*>(node);

                for (size_t i = 0; i < leaf->getLeafSize(); i++)
                {
                    update_rings(leaf->getDataset()->getElement(i), node);
                }

            }

        }

        void merge_rings_and_covarage(Node<O, T>* node)
        {

            if (node == nullptr)
                return;

            update_rings(node);
            Node<O, T>* update_ptr = node->getParent();

            while (update_ptr != nullptr)
            {
                update_covarage(node);
                update_rings(update_ptr);
                update_ptr = update_ptr->getParent();
            }

        }

        void update_covarage(Node<O, T>* node)
        {
            for (size_t i = 0; i < node->getChildrensSize(); i++)
            {
                node->setCoverage(std::max(node->getCoverage(), node->getChild(i)->getCoverage()));
            }
        }

        void update_rings(dataset::BasicArrayObject<O, T>& element, Node<O, T>* node)
        {

            double dist;

            for (size_t i = 0; i < this->numPivots; i++)
            {
                dist = this->distanceFunction->getDistance(element, this->pivots->getPivot(i));
                node->setRingMin(i, std::min(dist, node->getRingMin(i)));
                node->setRingMax(i, std::max(dist, node->getRingMax(i)));
            }

        }

        void split(dataset::BasicArrayObject<O, T>& element, Node<O, T>* node)
        {

            if (node->isLeafNode())
            {

                auto* leafNode = dynamic_cast<LeafNode<O, T>*>(node);

                leafNode->insert(element);
                std::unique_ptr<pivots::Pivot<O, T>> newPivots = pivots::PivotFactory<O, T>::createPivot(this->pivots->getPivotType(), this->pivots);
                newPivots->operator()(leafNode->getDataset(), this->distanceFunction, 2);

                std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>();
                std::unique_ptr<dataset::Dataset<O, T>> rightDataset = std::make_unique<dataset::Dataset<O, T>>();

                double dist1, dist2, cov1 = std::numeric_limits<double>::min(), cov2 = std::numeric_limits<double>::min();

                for (size_t i = 0; i < leafNode->getLeafSize(); i++)
                {

                    dist1 = this->distanceFunction->operator()(newPivots->getPivot(0), leafNode->getDataset()->getElement(i));
                    dist2 = this->distanceFunction->operator()(newPivots->getPivot(1), leafNode->getDataset()->getElement(i));

                    if (dist1 < dist2) {
                        cov1 = std::max(cov1, dist1);
                        leftDataset->insert(leafNode->getDataset()->getElement(i));
                    }
                    else {
                        cov2 = std::max(cov2, dist2);
                        rightDataset->insert(leafNode->getDataset()->getElement(i));
                    }
                }

                if (node->getParent() == nullptr) //root
                {
                    std::unique_ptr<DirectoryNode<O, T>> newDirectoryNode = std::make_unique<DirectoryNode<O, T>>(this->numPivots);
                    newDirectoryNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(0)));

                    std::unique_ptr<LeafNode<O, T>> newLeafNode2 = std::make_unique<LeafNode<O, T>>(this->numPivots, std::move(rightDataset));
                    newLeafNode2->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(1)));

                    leafNode->getDataset()->clear();
                    leafNode->getDataset().reset();
                    leafNode->setDataset(std::move(leftDataset));
                    leafNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(0)));

                    newDirectoryNode->addChild(std::move(root));
                    newDirectoryNode->addChild(std::move(newLeafNode2));

                    newDirectoryNode->getChild(1)->setCoverage(cov1);
                    newDirectoryNode->getChild(1)->setCoverage(cov2);

                    root = std::move(newDirectoryNode);

                    double dist_p1 = this->distanceFunction->getDistance(*root->getChild(0)->getPivot(), *root->getPivot()),
                           dist_p2 = this->distanceFunction->getDistance(*root->getChild(1)->getPivot(), *root->getPivot());

                    root->getChild(0)->setDistanceToParent(dist_p1);
                    root->getChild(1)->setDistanceToParent(dist_p2);
                    root->setCoverage(std::max(cov1 + dist_p1, cov2 + dist_p2));

                    update_rings(root->getChild(0).get());
                    update_rings(root->getChild(1).get());
                    merge_rings_and_covarage(root->getChild(0).get());

                }
                else
                {
                    std::unique_ptr<LeafNode<O, T>> newLeafNode = std::make_unique<LeafNode<O, T>>(this->numPivots, std::move(rightDataset));
                    newLeafNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(1)));

                    leafNode->getDataset()->clear();
                    leafNode->getDataset().reset();
                    leafNode->setDataset(std::move(leftDataset));
                    leafNode->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(0)));

                    Node<O, T>* parent = leafNode->getParent();
                    leafNode->setDistanceToParent(this->distanceFunction->operator()(*leafNode->getPivot(), *parent->getPivot()));
                    leafNode->setCoverage(cov1);

                    newLeafNode->setDistanceToParent(this->distanceFunction->operator()(*newLeafNode->getPivot(), *parent->getPivot()));
                    newLeafNode->setCoverage(cov2);

                    parent->setCoverage(std::max(cov1 + leafNode->getDistanceToParent(), cov2 + newLeafNode->getDistanceToParent()));

                    if ((parent->getChildrensSize() + 1) > this->numPerLeaf)
                    {
                        parent->addChild(std::move(newLeafNode));
                        split(element, parent);
                    }
                    else
                    {
                        update_rings(newLeafNode.get());
                        parent->addChild(std::move(newLeafNode));
                        update_rings(leafNode);
                        merge_rings_and_covarage(leafNode);
                    }

                }

                newPivots->clear();
                newPivots.reset();

            }
            else
            {
                std::unique_ptr<dataset::Dataset<O, T>> pivotsData = std::make_unique<dataset::Dataset<O, T>>();
                std::unique_ptr<pivots::Pivot<O, T>> newPivots = pivots::PivotFactory<O, T>::createPivot(this->pivots->getPivotType(), this->pivots);

                for (size_t i = 0; i < node->getChildrensSize(); i++)
                {
                    pivotsData->insert(*node->getChild(i)->getPivot());
                }

                newPivots->operator()(pivotsData, this->distanceFunction, 2);

                double dist1, dist2, cov1, cov2;
                size_t pos_P1 = std::numeric_limits<size_t>::max(), pos_P2 = std::numeric_limits<size_t>::max();
                std::vector<std::unique_ptr<Node<O, T>>> childs = std::move(node->getChildrens()), left, right;

                for (size_t i = 0; i < childs.size(); i++)
                {
                    dist1 = this->distanceFunction->operator()(newPivots->getPivot(0), *childs.at(i)->getPivot());
                    dist2 = this->distanceFunction->operator()(newPivots->getPivot(1), *childs.at(i)->getPivot());

                    if (dist1 == 0.0)
                        pos_P1 = left.size();

                    if (dist2 == 0.0)
                        pos_P2 = right.size();

                    if (dist1 < dist2) {
                        cov1 = std::max(cov1, dist1 + childs.at(i)->getCoverage());
                        left.push_back(std::move(childs.at(i)));
                    }
                    else {
                        cov2 = std::max(cov2, dist2 + childs.at(i)->getCoverage());
                        right.push_back(std::move(childs.at(i)));
                    }
                }

                if (pos_P1 == std::numeric_limits<size_t>::max() || pos_P2 == std::numeric_limits<size_t>::max())
                    throw std::runtime_error("Error in split function");

                if (node->getParent() == nullptr) //root
                {
                    std::unique_ptr<Node<O, T>> newDirectoryNode = std::move(left[pos_P1]), newDirectoryNode2 = std::move(right[pos_P2]);

                    for (size_t i = 0; i < left.size(); i++)
                    {
                        if (i != pos_P1)
                            newDirectoryNode->addChild(std::move(left[i]));
                    }

                    for (size_t i = 0; i < right.size(); i++)
                    {
                        if (i != pos_P2)
                            newDirectoryNode2->addChild(std::move(right[i]));
                    }

                    newDirectoryNode->setCoverage(cov1);
                    newDirectoryNode2->setCoverage(cov2);

                    update_rings(newDirectoryNode.get());
                    update_rings(newDirectoryNode2.get());

                    root->clear();
                    root.reset();
                    root = std::make_unique<DirectoryNode<O, T>>(this->numPivots);
                    root->setPivot(std::make_unique<dataset::BasicArrayObject<O, T>>(newPivots->getPivot(0)));
                    root->addChild(std::move(newDirectoryNode));
                    root->addChild(std::move(newDirectoryNode2));

                    double dist_p1 = this->distanceFunction->getDistance(*root->getChild(0)->getPivot(), *root->getPivot()),
                           dist_p2 = this->distanceFunction->getDistance(*root->getChild(1)->getPivot(), *root->getPivot());

                    root->getChild(0)->setDistanceToParent(dist_p1);
                    root->getChild(1)->setDistanceToParent(dist_p2);

                    root->setCoverage(std::max(cov1 + dist_p1, cov2 + dist_p2));
                    update_rings(root.get());

                }
                else
                {
                    throw std::runtime_error("Not implemented yet");
                }

            }

        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("Not implemented yet");
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("Not implemented yet");
        }

    };

}

#endif //GERVLIB_PMTREE_H
