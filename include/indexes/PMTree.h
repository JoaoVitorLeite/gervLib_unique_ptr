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
            memoryStatus = MEMORY_STATUS::NONE;
        }

        explicit Node(size_t num_pivots)
        {
            parent = nullptr;
            pivot = nullptr;
            dist_to_parent = 0.0;
            coverage = 0.0;
            nodeID = 0;
            rings = std::make_unique<std::vector<std::pair<double, double>>>(num_pivots, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
            memoryStatus = MEMORY_STATUS::NONE;
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
            utils::check_range(0, childrens->size()-1, index, "Index out of bounds");
            childrens[index] = std::move(child);
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

}

#endif //GERVLIB_PMTREE_H
