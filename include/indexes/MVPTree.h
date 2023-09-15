//
// Created by joaovictor on 25/08/23.
//

#ifndef GERVLIB_MVPTREE_H
#define GERVLIB_MVPTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "NAryTree.h"
#include <queue>
#include <map>

namespace gervLib::index::mvptree
{

    template <typename O, typename T>
    class Node : public serialize::Serialize
    {

    protected:
        std::unique_ptr<std::vector<dataset::BasicArrayObject<O, T>>> pivots;
        std::vector<std::unique_ptr<Node<O, T>>> childrens;
        size_t nodeID{};
        std::unique_ptr<std::vector<std::vector<double>>> splits;
        MEMORY_STATUS memoryStatus = MEMORY_STATUS::NONE;

    public:
        Node() = default;

        explicit Node(size_t fo)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
        }

        Node(size_t lpn, size_t fo)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
            pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(lpn);
        }

        Node(size_t lpn, size_t fo, size_t ns)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
            pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(lpn);
            splits = std::make_unique<std::vector<std::vector<double>>>(lpn, std::vector<double>(ns, -1.0));
        }

        virtual ~Node()
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }

            if (splits != nullptr)
            {
                splits->clear();
                splits.reset();
            }
        }

        void setNumberOfPivots(size_t lpn)
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }
            pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(lpn);
        }

        void setNumberOfFanout(size_t fo)
        {
            if (!childrens.empty())
                childrens.clear();
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
        }

        void setNumberOfSplits(size_t lpn, size_t ns)
        {
            if (splits != nullptr)
            {
                splits->clear();
                splits.reset();
            }
            splits = std::make_unique<std::vector<std::vector<double>>>(lpn, std::vector<double>(ns, -1.0));
        }

        void setChildren(size_t index, std::unique_ptr<Node<O, T>> node)
        {
            utils::check_range(0, childrens.size() - 1, index, "Node::setChildren(): Index out of range.");
            childrens[index] = std::move(node);
        }

        void setPivot(size_t index, dataset::BasicArrayObject<O, T> pivot)
        {

            if (pivots == nullptr)
                throw std::runtime_error("Node::setPivot(): Pivots vector is null.");

            utils::check_range(0, childrens.size() - 1, index, "Node::setChildren(): Index out of range.");
            pivots->at(index) = pivot;
        }

        void setSplit(size_t r, size_t c, double value)
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::setSplit(): Splits vector is null.");

            utils::check_range(0, splits->size() - 1, r, "Node::setSplit(): Index out of range.");
            utils::check_range(0, splits->at(r).size() - 1, c, "Node::setSplit(): Index out of range.");

            splits->at(r)[c] = value;

        }

        void setNodeID(size_t id) { nodeID = id; }

        void setMemoryStatus(MEMORY_STATUS ms) { memoryStatus = ms; }

        std::unique_ptr<Node<O, T>>& getChildren(size_t index)
        {
            utils::check_range(0, childrens.size() - 1, index, "Node::getChildren(): Index out of range.");
            return childrens[index];
        }

        dataset::BasicArrayObject<O, T>& getPivot(size_t index)
        {
            if (pivots == nullptr)
                throw std::runtime_error("Node::getPivot(): Pivots vector is null.");

            utils::check_range(0, pivots->size() - 1, index, "Node::getPivot(): Index out of range.");
            return pivots->at(index);
        }

        double getSplit(size_t r, size_t c)
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::getSplit(): Splits vector is null.");

            utils::check_range(0, splits->size() - 1, r, "Node::getSplit(): Index out of range.");
            utils::check_range(0, splits->at(r).size() - 1, c, "Node::getSplit(): Index out of range.");

            return splits->at(r)[c];
        }

        MEMORY_STATUS getMemoryStatus() { return memoryStatus; }

        size_t getNumberOfPivots()
        {
            if (pivots == nullptr)
                throw std::runtime_error("Node::getNumberOfPivots(): Pivots vector is null.");

            return pivots->size();
        }

        size_t getNumberOfSplits()
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::getNumberOfSplits(): Splits vector is null.");

            return splits->at(0).size();
        }

        size_t getNodeID() { return nodeID; }

        bool operator==(Node<O, T> *other)
        {
            return isEqual(other);
        }

        bool operator!=(Node<O, T> *other)
        {
            return !isEqual(other);
        }

        //virtual methods

        virtual void clear()
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }

            if (splits != nullptr)
            {
                splits->clear();
                splits.reset();
            }
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

                if (((pivots == nullptr && other->pivots != nullptr) || (pivots != nullptr && other->pivots == nullptr)) || (pivots->size() != other->pivots->size()))
                    return false;
                else
                {
                    for (size_t i = 0; i < pivots->size(); i++)
                    {
                        if (pivots->at(i) != other->pivots->at(i))
                            return false;
                    }
                }

                if (((splits == nullptr && other->splits != nullptr) || (splits != nullptr && other->splits == nullptr)) || (splits->size() != other->splits->size() || splits->at(0).size() != other->splits->at(0).size()))
                    return false;
                else
                {
                    for (size_t i = 0; i < splits->size(); i++)
                    {
                        if (splits->at(i) != other->splits->at(i))
                            return false;
                    }
                }

                return true;

            }
        }

        virtual void print(std::ostream& os) const
        {
            os << "Node type: Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;

            if (pivots != nullptr)
            {
                os << "Pivots: " << std::endl;
                for (size_t i = 0; i < pivots->size(); i++)
                    os << pivots->at(i) << std::endl;
            }
            else
                os << "Pivots: NULL" << std::endl;

            if (splits != nullptr)
            {
                os << "Splits: " << std::endl;
                for (size_t i = 0; i < splits->size(); i++)
                {
                    for (size_t j = 0; j < splits->at(i).size(); j++)
                        os << splits->at(i)[j] << " ";
                    os << std::endl;
                }
            }
            else
                os << "Splits: NULL" << std::endl;

            os << "Memory Status: " << gervLib::index::memoryStatusMap[memoryStatus] << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            if (pivots != nullptr)
            {
                sz = pivots->size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                for (auto pvt : *pivots)
                {
                    sz = pvt.getSerializedSize();
                    memcpy(data.get() + offset, &sz, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> pvtData = pvt.serialize();
                    memcpy(data.get() + offset, pvtData.get(), sz);
                    offset += sz;
                    pvtData.reset();
                }
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (splits != nullptr)
            {
                sz = splits->size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                sz = splits->at(0).size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                for (auto split : *splits)
                {
                    memcpy(data.get() + offset, split.data(), sizeof(double) * split.size());
                    offset += sizeof(double) * split.size();
                }
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            std::string aux = gervLib::index::memoryStatusMap[memoryStatus];
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

            memcpy(&nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            size_t sz_loop;
            memcpy(&sz_loop, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz_loop != 0)
            {

                if (pivots != nullptr)
                {
                    pivots->clear();
                    pivots.reset();
                }

                pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(sz_loop);

                for (size_t i = 0; i < sz_loop; i++)
                {
                    memcpy(&sz, _data.get() + offset, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> pvtData(new u_char[sz]);
                    memcpy(pvtData.get(), _data.get() + offset, sz);
                    offset += sz;

                    pivots->at(i).deserialize(std::move(pvtData));
                }
            }
            else
                pivots = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            size_t sz2;
            memcpy(&sz2, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0 && sz2 != 0)
            {

                if (splits != nullptr)
                {
                    splits->clear();
                    splits.reset();
                }

                splits = std::make_unique<std::vector<std::vector<double>>>(sz, std::vector<double>(sz2));

                for (size_t i = 0; i < sz; i++)
                {
                    memcpy(splits->at(i).data(), _data.get() + offset, sizeof(double) * sz2);
                    offset += sizeof(double) * sz2;
                }
            }
            else
                splits = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(aux.data(), _data.get() + offset, sz);
            offset += sz;

            memoryStatus = gervLib::index::memoryStatusMapReverse[aux];

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            ans += sizeof(size_t); //nodeID

            ans += sizeof(size_t); //pivots->size()
            if (pivots != nullptr)
            {

                for (auto pvt : *pivots)
                    ans += sizeof(size_t) + pvt.getSerializedSize();

            }

            ans += sizeof(size_t); //splits->size()
            ans += sizeof(size_t); //splits->at(0).size()
            if (splits != nullptr)
            {
                ans += sizeof(double) * splits->size() * splits->at(0).size();
            }

            ans += sizeof(size_t) + memoryStatusMap[memoryStatus].size(); //memoryStatus

            return ans;

        }

    };

    template <typename O, typename T>
    class DirectoryNode : public Node<O, T>
    {

    public:
        DirectoryNode() = default;

        explicit DirectoryNode(size_t fo) : Node<O, T>(fo) {}

        DirectoryNode(size_t lpn, size_t fo) : Node<O, T>(lpn, fo) {}

        DirectoryNode(size_t lpn, size_t fo, size_t ns) : Node<O, T>(lpn, fo, ns) {}

        ~DirectoryNode() override = default;

        bool isLeafNode() override
        {
            return false;
        }

        void print(std::ostream& os) const override
        {
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;

            if (this->pivots != nullptr)
            {
                os << "Pivots: " << std::endl;
                for (size_t i = 0; i < this->pivots->size(); i++)
                    os << this->pivots->at(i) << std::endl;
            }
            else
                os << "Pivots: NULL" << std::endl;

            if (this->splits != nullptr)
            {
                os << "Splits: " << std::endl;
                for (size_t i = 0; i < this->splits->size(); i++)
                {
                    for (size_t j = 0; j < this->splits->at(i).size(); j++)
                        os << this->splits->at(i)[j] << " ";
                    os << std::endl;
                }
            }
            else
                os << "Splits: NULL" << std::endl;

            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;

        }

    };

    template <typename O, typename T>
    class LeafNode : public Node<O, T> {

    private:
        std::unique_ptr<dataset::Dataset<O, T>> dataset;
        std::unique_ptr<Index<O, T>> index;

    public:
        LeafNode() : Node<O, T>() {
            this->dataset = nullptr;
            this->index = nullptr;
        }

        explicit LeafNode(size_t lpn) : Node<O, T>(lpn, 0) {
            this->dataset = nullptr;
            this->index = nullptr;
        }

        explicit LeafNode(size_t lpn, std::unique_ptr<dataset::Dataset<O, T>> _dataset) : Node<O, T>(lpn, 0) {
            this->dataset = std::move(_dataset);
            this->index = nullptr;
        }

        LeafNode(size_t lpn, std::unique_ptr<dataset::Dataset<O, T>> _dataset, std::unique_ptr<Index<O, T>> _index) : Node<O, T>(lpn, 0) {
            this->dataset = std::move(_dataset);
            this->index = std::move(_index);
        }

        ~LeafNode() override {
            Node<O, T>::clear();
            if (this->dataset != nullptr)
                this->dataset.reset();
            if (this->index != nullptr)
                this->index.reset();
        }

        bool isLeafNode() override {
            return true;
        }

        void clear() override {
            Node<O, T>::clear();
            if (this->dataset != nullptr)
                this->dataset.reset();
            if (this->index != nullptr)
                this->index.reset();
        }

        void setIndex(std::unique_ptr<Index<O, T>> _index) {
            this->index.reset();
            this->index = std::move(_index);
        }

        void setDataset(std::unique_ptr<dataset::Dataset<O, T>> _dataset) {
            this->dataset.reset();
            this->dataset = std::move(_dataset);
        }

        std::unique_ptr<Index<O, T>> &getIndex() {
            return this->index;
        }

        std::unique_ptr<dataset::Dataset<O, T>> &getDataset() {
            return this->dataset;
        }

        void insert(std::unique_ptr<dataset::BasicArrayObject<O, T>> &object) {
            if (dataset == nullptr)
                dataset = std::make_unique<dataset::Dataset<O, T>>();
            this->dataset->insert(*object);
        }

        bool isEqual(std::unique_ptr<Node<O, T>>& other) override {

            auto* node = dynamic_cast<LeafNode<O, T>*>(other.get());

            if (node == nullptr)
                return false;

            if (this->memoryStatus == MEMORY_STATUS::IN_DISK || node->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                return this->nodeID == node->getNodeID();
            else {

                if (((this->pivots == nullptr && node->getNumberOfPivots() != 0) || (this->pivots != nullptr && node->pivots == nullptr)) || (this->pivots->size() != node->pivots->size()))
                    return false;
                else {
                    for (size_t i = 0; i < this->pivots->size(); i++) {
                        if (this->pivots->at(i) != node->pivots->at(i))
                            return false;
                    }
                }

                if (((this->splits == nullptr && node->splits != nullptr) || (this->splits != nullptr && node->splits == nullptr)) || (this->splits->size() != node->splits->size() || this->splits->at(0).size() != node->splits->at(0).size()))
                    return false;
                else {
                    for (size_t i = 0; i < this->splits->size(); i++) {
                        if (this->splits->at(i) != node->splits->at(i))
                            return false;
                    }
                }

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

                return true;

            }

        }

        void print(std::ostream& os) const override
        {
            os << "Node type: Leaf Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;

            if (this->pivots != nullptr)
            {
                os << "Pivots: " << std::endl;
                for (size_t i = 0; i < this->pivots->size(); i++)
                    os << this->pivots->at(i) << std::endl;
            }
            else
                os << "Pivots: NULL" << std::endl;

            if (this->splits != nullptr)
            {
                os << "Splits: " << std::endl;
                for (size_t i = 0; i < this->splits->size(); i++)
                {
                    for (size_t j = 0; j < this->splits->at(i).size(); j++)
                        os << this->splits->at(i)[j] << " ";
                    os << std::endl;
                }
            }
            else
                os << "Splits: NULL" << std::endl;

            os << "Memory Status: " << gervLib::index::memoryStatusMap[this->memoryStatus] << std::endl;

            if (this->dataset != nullptr)
                os << "Dataset: " << std::endl << *this->dataset << std::endl;
            else
                os << "Dataset: NULL" << std::endl;

            if (this->index != nullptr)
                os << "Index: " << std::endl << *this->index << std::endl;
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
    class MVPTree : public Index<O, T>
    {

    private:
        std::unique_ptr<Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{}, branchFactor{}, levelsPerNode{}, fanout{}, numMaxSplits{};
        bool storePivotsInLeaf{}, storeDirectoryNode{}, storeLeafNode{}, useLAESA{};
        std::string serializedTree{};

    private:
        enum cases{R0, R1, R2, R3, R4};

    private:
        std::vector<double> calculatePivotDistances(dataset::BasicArrayObject<O, T>& vp, std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {
            std::vector<double> ans = std::vector<double>(dataset->getCardinality(), 0.0);

            for (size_t i = 0; i < dataset->size(); i++)
                ans[i] = this->distanceFunction->operator()(vp, dataset->getElement(i));

            return ans;
        }

        void calculateSplits(Node<O, T>* node, std::vector<double>& dists, size_t n, size_t split_index)
        {
            size_t len = branchFactor - 1;

            if (!dists.empty())
            {
                if (node->getSplit(n, split_index * len) == -1.0)
                {
                    std::vector<double> tmp = dists;
                    std::sort(tmp.begin(), tmp.end());
                    double factor = (tmp.size() * 1.0)/branchFactor;

                    for(size_t i = 0; i < len; i++)
                    {

                        size_t pos = (i + 1) * factor;
                        size_t lo = floor(pos);
                        size_t hi = (pos <= tmp.size()) ? ceil(pos) : 0;
                        node->setSplit(n, split_index * len + i, (tmp[lo] + tmp[hi])/2.0);

                    }

                }

            }

        }

        std::unique_ptr<dataset::Dataset<O, T>> cullPoints(std::unique_ptr<dataset::Dataset<O, T>>& dataset, std::vector<double> dists, double m, bool less)
        {
            std::unique_ptr<dataset::Dataset<O, T>> ans = std::make_unique<dataset::Dataset<O, T>>();
            ans->setPath(dataset->getPath());
            ans->setSeed(dataset->getSeed());
            ans->setDimensionality(dataset->getDimensionality());

            for(size_t i = 0; i < dataset->getCardinality(); i++)
            {

                if(less)
                {

                    if(dists[i] <= m)
                    {

                        ans->insert(dataset->getElement(i));

                    }

                }
                else
                {

                    if(dists[i] > m)
                    {

                        ans->insert(dataset->getElement(i));

                    }

                }

            }

            if (ans->getCardinality() == 0)
                ans.reset();

            return ans;
        }

        std::vector<std::unique_ptr<dataset::Dataset<O, T>>> splitNode(Node<O, T>* node, std::unique_ptr<dataset::Dataset<O, T>> dataset)
        {

            std::map<long long, std::unique_ptr<dataset::Dataset<O, T>>> pnts, pnts2;
            pnts[0] = std::move(dataset);

            std::vector<std::unique_ptr<dataset::Dataset<O, T>>> ans(this->fanout);

            size_t len = branchFactor - 1, n = 0;
            double m;
            std::unique_ptr<dataset::Dataset<O, T>> culledPts = nullptr;

            do{

                for(auto it = pnts.begin(); it != pnts.end(); it++)
                {

                    long long node_index = it->first;
                    std::unique_ptr<dataset::Dataset<O, T>> node_dataset = std::move(it->second);
                    dataset::BasicArrayObject<O, T> vp = node->getPivot(n);
                    std::vector<double> dists = calculatePivotDistances(vp, node_dataset);

                    if (!dists.empty())
                    {
                        calculateSplits(node, dists, n, node_index);

                        for(size_t j = 0; j < len; j++)
                        {
                            m = node->getSplit(n, node_index * len + j);
                            culledPts = cullPoints(node_dataset, dists, m, true);

                            if(culledPts != nullptr)
                            {

                                pnts2[node_index * branchFactor + j] = std::move(culledPts);

                            }

                        }

                        m = node->getSplit(n, node_index * len + len - 1);
                        culledPts = cullPoints(node_dataset, dists, m, false);

                        if(culledPts != nullptr)
                        {

                            pnts2[node_index * branchFactor + branchFactor - 1] = std::move(culledPts);

                        }

                    }

                    node_dataset->clear();
                    node_dataset.reset();

                }

                pnts = std::move(pnts2);
                n++;

            } while (n < levelsPerNode);

            for(auto it = pnts.begin(); it != pnts.end(); it++)
            {

                long long index = it->first;
                std::unique_ptr<dataset::Dataset<O, T>> list = std::move(it->second);

                if(list != nullptr)
                {

                    ans[index] = std::move(list);

                }

            }

            return ans;

        }

        void deleteRecursive(std::unique_ptr<Node<O, T>> node)
        {

            if (node == nullptr)
                return;

            if (!node->isLeafNode()) {

                for (size_t i = 0; i < fanout; i++)
                    deleteRecursive(std::move(node->getChildren(i)));

            }

            node->clear();
            node.reset();

        }

        void clearRecursive(std::unique_ptr<Node<O, T>>& node)
        {

            if (node == nullptr)
                return;

            if (!node->isLeafNode()) {

                for (size_t i = 0; i < fanout; i++)
                    clearRecursive(node->getChildren(i));

            }

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

            if (!node1->isLeafNode())
            {
                for (size_t i = 0; i < fanout; i++)
                {
                    if (!isEqualHelper(node1->getChildren(i), node2->getChildren(i)))
                        return false;
                }
            }

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

            std::string result = (node->isLeafNode() ? "L" : "D") + std::to_string(node->getNodeID()) + " " + std::to_string(fanout) + " ";

            if (node->isLeafNode())
                for (size_t i = 0; i < fanout; i++)
                    result += "null ";
            else
                for (size_t i = 0; i < fanout; i++)
                    result += serializeTreeRecursive(node->getChildren(i));

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

                node->setNumberOfFanout(fanout);

                for (size_t i = 0; i < fanout; i++)
                {
                    buildTree(node->getChildren(i), aux->children[i]);
                }

            }

        }

        double minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3)
        {
            double ans{};

            if (sqCase == R0)
            {

                if(nodeCase == R0)
                    ans = 0.0;
                else if(nodeCase == R1)
                    ans = fabs(d_sq_p2 - mu2);
                else if(nodeCase == R2)
                {

                    if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3)))
                    {
                        ans = fabs(d_sq_p2 - mu3);
                    }
                    else
                    {
                        ans = fabs(d_sq_p1 - mu1);

                    }

                }
                else if(nodeCase == R3)
                {

                    if ((d_sq_p2 + d_sq_p1 + mu1) <= mu3)
                    {
                        ans = fabs(d_sq_p2 - mu3);
                    }
                    else
                    {
                        if ((d_sq_p2 + d_sq_p1 + mu2) <= mu3)
                        {
                            ans = fabs(d_sq_p2 - mu3);
                        } else {
                            ans = fabs(d_sq_p1 - mu1);
                        }
                    }

                }
                else
                    throw std::runtime_error("non-existent case");

            }
            else if(sqCase == R1)
            {

                if(nodeCase == R0)
                    ans = fabs(d_sq_p2 - mu2);
                else if(nodeCase == R1)
                    ans = 0.0;
                else if(nodeCase == R2)
                {

                    if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3)))
                    {
                        ans = fabs(d_sq_p2 - mu3);
                    }
                    else
                    {
                        ans = fabs(d_sq_p1 - mu1);
                    }

                }
                else if(nodeCase == R3)
                {

                    if ((d_sq_p2 + d_sq_p1 + mu1) <= mu3)
                    {
                        ans = fabs(d_sq_p2 - mu3);
                    }
                    else
                    {
                        ans = fabs(d_sq_p1 - mu1);

                    }

                }
                else
                    throw std::runtime_error("non-existent case");

            }
            else if(sqCase == R2)
            {

                if(nodeCase == R0)
                {
                    if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){
                        ans = fabs(d_sq_p2 - mu2);
                    } else {
                        ans = fabs(d_sq_p1 - mu1);
                    }
                }
                else if(nodeCase == R1)
                {
                    ans = fabs(d_sq_p1 - mu1);
                }
                else if(nodeCase == R2)
                    ans = 0.0;
                else if(nodeCase == R3)
                {
                    ans = fabs(d_sq_p2 - mu3);
                }
                else
                    throw std::runtime_error("non-existent case");

            }
            else if(sqCase == R3)
            {

                if(nodeCase == R0)
                {
                    if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){
                        ans = fabs(d_sq_p2 - mu2);
                    } else {
                        ans = fabs(d_sq_p1 - mu1);
                    }

                }
                else if(nodeCase == R1){
                    ans = fabs(d_sq_p1 - mu1);
                }
                else if(nodeCase == R2){
                    ans = fabs(d_sq_p2 - mu3);
                }
                else if(nodeCase == R3)
                    ans = 0.0;
                else
                    throw std::runtime_error("non-existent case");

            }
            else if(sqCase == R4)
            {

                if(nodeCase == R0)
                {
                    if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){ //mu2 está contido em mu1
                        ans = fabs(d_sq_p2 - mu2);
                    } else {
                        ans = fabs(d_sq_p1 - mu1);
                    }
                }
                else if(nodeCase == R1){
                    ans = fabs(d_sq_p1 - mu1);
                }
                else if(nodeCase == R2){
                    ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
                }
                else if(nodeCase == R3){
                    ans = fabs(d_sq_p1 - M);

                }
                else
                    throw std::runtime_error("non-existent case");

            }

            return ans;

        }

        double maxDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3)
        {

            double ans{};

            if (sqCase == R0)
            {
                if (nodeCase == R0)
                    ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
                else if (nodeCase == R1)
                    ans = d_sq_p1 + mu1;
                else if (nodeCase == R2)
                    ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
                else if (nodeCase == R3)
                    ans = d_sq_p1 + M;
                else
                    throw std::runtime_error("non-existent case");
            }
            else if (sqCase == R1)
            {
                if (nodeCase == R0)
                    ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
                else if (nodeCase == R1)
                    ans = d_sq_p1 + mu1;
                else if (nodeCase == R2)
                    ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
                else if (nodeCase == R3)
                    ans = d_sq_p1 + M;
                else
                    throw std::runtime_error("non-existent case");
            }
            else if (sqCase == R2)
            {
                if (nodeCase == R0)
                    ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
                else if (nodeCase == R1)
                    ans = d_sq_p1 + mu1;
                else if (nodeCase == R2)
                    ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
                else if (nodeCase == R3)
                    ans = d_sq_p1 + M;
                else
                    throw std::runtime_error("non-existent case");
            }
            else if (sqCase == R3)
            {
                if (nodeCase == R0)
                    ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
                else if (nodeCase == R1)
                    ans = d_sq_p1 + mu1;
                else if (nodeCase == R2)
                    ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
                else if (nodeCase == R3)
                    ans = d_sq_p1 + M;
                else
                    throw std::runtime_error("non-existent case");
            }
            else if (sqCase == R4)
            {
                if (nodeCase == R0)
                    ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
                else if (nodeCase == R1)
                    ans = d_sq_p1 + mu1;
                else if (nodeCase == R2)
                    ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
                else if (nodeCase == R3)
                    ans = d_sq_p1 + M;
                else
                    throw std::runtime_error("non-existent case");
            }

            return ans;

        }

    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount,iowrite,ioread,seed";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning,iowrite,ioread";
        }

    public:
        MVPTree()
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
            this->branchFactor = 0;
            this->levelsPerNode = 0;
            this->fanout = 0;
            this->numMaxSplits = 0;
            this->storePivotsInLeaf = false;
            this->storeDirectoryNode = false;
            this->storeLeafNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::MVPTREE;
            this->indexName = "MVPTREE";
            this->indexFolder = "";
        }

        MVPTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
                std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
                std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0, size_t _branchFactor = 2,
                size_t _levelsPerNode = 2, size_t _fanout = 4, size_t _numMaxSplits = 2, bool _storePivotsInLeaf = true, bool _storeDirectoryNode = false,
                bool _storeLeafNode = false, bool _useLAESA = true, std::string folder = "")
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
            this->branchFactor = _branchFactor;
            this->levelsPerNode = _levelsPerNode;
            this->fanout = _fanout;
            this->numMaxSplits = _numMaxSplits;
            this->storePivotsInLeaf = _storePivotsInLeaf;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->storeLeafNode = _storeLeafNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::MVPTREE;
            this->indexName = "MVPTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("mvp_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit MVPTree(std::string _folder, std::string serializedFile = "")
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
            this->branchFactor = 0;
            this->levelsPerNode = 0;
            this->fanout = 0;
            this->numMaxSplits = 0;
            this->storePivotsInLeaf = false;
            this->storeDirectoryNode = false;
            this->storeLeafNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::MVPTREE;
            this->indexName = "MVPTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);
        }

        ~MVPTree()
        {
            deleteRecursive(std::move(root));
        }

        void clear() override
        {
            clearRecursive(root);
        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<MVPTree<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

        }

        std::unique_ptr<Node<O, T>>& getRoot()
        {
            return this->root;
        }

        void print(std::ostream& os) const override {

            std::stack<std::pair<Node < O, T> *, size_t>>
            nodeStack;
            nodeStack.push(std::make_pair(root.get(), 0));

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "MVPTree" << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Number of pivots per node: " << levelsPerNode << std::endl;
            os << "Number of points per leaf: " << numPerLeaf << std::endl;
            os << "Branching factor: " << branchFactor << std::endl;
            os << "Fanout: " << fanout << std::endl;
            os << "Number of max splits: " << numMaxSplits << std::endl;
            os << "Store pivots in leaf: " << (storePivotsInLeaf ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Use LAESA: " << (useLAESA ? "true" : "false") << std::endl;

            while (!nodeStack.empty())
            {

                auto currentNode = nodeStack.top();
                nodeStack.pop();

                os << "**********************************************************************************************************************************************************************************\n\n";
                os << *currentNode.first << std::endl;

                if (!currentNode.first->isLeafNode()) {

                    for (size_t i = 0; i < fanout; i++) {
                        if (currentNode.first->getChildren(i) != nullptr)
                            nodeStack.push(
                                    std::make_pair(currentNode.first->getChildren(i).get(), currentNode.second + 1));
                    }

                }

            }

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
            std::unique_ptr<pivots::Pivot<O, T>> globalPivots;

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

                currentNode.first->setNumberOfPivots(levelsPerNode);
                this->pivots->operator()(currentNode.second, this->distanceFunction, levelsPerNode);

                for (size_t i = 0; i < levelsPerNode; i++) {
                    currentNode.first->setPivot(i, this->pivots->getPivot(i));

                    if (!storePivotsInLeaf)
                    {
                        currentNode.second->erase(this->pivots->getPivot(i));
                    }

                }

                if (currentNode.second->getCardinality() <= (numPerLeaf + levelsPerNode))
                {
                    auto* leafNode = (LeafNode<O, T>*) currentNode.first;
                    leafNode->setNodeID(currentNodeID++);

                    if (useLAESA) {
                        std::filesystem::path leafIndexPath(this->indexFolder);
                        leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                this->distanceFunction->getDistanceType());
                        std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                std::move(currentNode.second), std::move(df),
                                pivots::PivotFactory<O, T>::clone(globalPivots), this->numPivots, leafIndexPath);
                        this->distanceFunction->setDistanceCount(this->distanceFunction->getDistanceCount() + idx->getDistanceCount());
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
                    currentNode.first->setNumberOfFanout(fanout);
                    currentNode.first->setNumberOfSplits(levelsPerNode, numMaxSplits);
                    std::vector<std::unique_ptr<dataset::Dataset<O, T>>> childrensDataset = splitNode(currentNode.first, std::move(currentNode.second));

                    for(size_t i = 0; i < fanout; i++)
                    {

                        if(childrensDataset[i] != nullptr)
                        {

                            if (childrensDataset[i]->getCardinality() <= (numPerLeaf + levelsPerNode))
                            {
                                currentNode.first->setChildren(i, std::make_unique<LeafNode<O, T>>());
                                nodeQueue.push(std::make_pair(currentNode.first->getChildren(i).get(), std::move(childrensDataset[i])));
                            }
                            else
                            {
                                currentNode.first->setChildren(i, std::make_unique<DirectoryNode<O, T>>());
                                nodeQueue.push(std::make_pair(currentNode.first->getChildren(i).get(), std::move(childrensDataset[i])));
                            }

                        }

                    }

                    auto *directoryNode = (DirectoryNode<O, T> *) currentNode.first;
                    directoryNode->setNodeID(currentNodeID++);

                    if (this->storeDirectoryNode) {
                        directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                        this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                        directoryNode->clear();
                    }
                    else
                        directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

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
                      << std::to_string(configure::IORead - ioR)
                      << ","
                      << std::to_string(this->pivots->getSeed()) << std::endl;
            buildFile.close();

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
            double dist, d_sq_p1, d_sq_p2;
            std::vector<cases> casesNodeVec = {R0, R1, R2, R3};

            nodeQueue.push(query::Partition<Node<O, T>*>(root.get(), 0.0, std::numeric_limits<double>::max()));

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty())) {

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

                            for (size_t i = 0 ; i < levelsPerNode; i++)
                            {
                                dist = this->distanceFunction->operator()(query, currentLeafNode->getPivot(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getPivot(i).getOID(), dist));
                            }

                        }

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

                    }
                    else
                    {

                        cases sqCase;
                        d_sq_p1 = this->distanceFunction->operator()(query, currentNode->getPivot(0));
                        d_sq_p2 = this->distanceFunction->operator()(query, currentNode->getPivot(1));

                        if (!storePivotsInLeaf)
                        {
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot(0).getOID(), d_sq_p1));
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot(1).getOID(), d_sq_p2));
                        }

                        if(d_sq_p1 <= currentNode->getSplit(0, 0) && d_sq_p2 <= currentNode->getSplit(1, 0))
                            sqCase = R0;
                        else if(d_sq_p1 <= currentNode->getSplit(0, 0) && d_sq_p2 > currentNode->getSplit(1, 0))
                            sqCase = R1;
                        else if(d_sq_p1 > currentNode->getSplit(0, 0) && d_sq_p2 <= currentNode->getSplit(1, 1))
                            sqCase = R2;
                        else if(d_sq_p1 > currentNode->getSplit(0, 0) && d_sq_p2 > currentNode->getSplit(1, 1))
                            sqCase = R3;
                        else
                            sqCase = R4;

                        for(size_t i = 0; i < fanout; i++)
                        {

                            if(currentNode->getChildren(i) != nullptr)
                            {

                                nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getChildren(i).get(),
                                                               minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, currentNode->getSplit(0, 0), currentPartition.getMax(), currentNode->getSplit(1, 0), currentNode->getSplit(1, 1)),
                                                               maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, currentNode->getSplit(0, 0), currentPartition.getMax(), currentNode->getSplit(1, 0), currentNode->getSplit(1, 1))));

                            }

                        }

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

                            for (size_t i = 0 ; i < levelsPerNode; i++)
                            {
                                dist = this->distanceFunction->operator()(query, currentLeafNode->getPivot(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getPivot(i).getOID(), dist));
                            }

                        }

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

                    }
                    else
                    {

                        cases sqCase;
                        d_sq_p1 = this->distanceFunction->operator()(query, currentNode->getPivot(0));
                        d_sq_p2 = this->distanceFunction->operator()(query, currentNode->getPivot(1));

                        if (!storePivotsInLeaf)
                        {
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot(0).getOID(), d_sq_p1));
                            elementQueue.push(query::ResultEntry<O>(currentNode->getPivot(1).getOID(), d_sq_p2));
                        }

                        if(d_sq_p1 <= currentNode->getSplit(0, 0) && d_sq_p2 <= currentNode->getSplit(1, 0))
                            sqCase = R0;
                        else if(d_sq_p1 <= currentNode->getSplit(0, 0) && d_sq_p2 > currentNode->getSplit(1, 0))
                            sqCase = R1;
                        else if(d_sq_p1 > currentNode->getSplit(0, 0) && d_sq_p2 <= currentNode->getSplit(1, 1))
                            sqCase = R2;
                        else if(d_sq_p1 > currentNode->getSplit(0, 0) && d_sq_p2 > currentNode->getSplit(1, 1))
                            sqCase = R3;
                        else
                            sqCase = R4;

                        for(size_t i = 0; i < fanout; i++)
                        {

                            if(currentNode->getChildren(i) != nullptr)
                            {

                                nodeQueue.push(query::Partition<Node<O, T>*>(currentNode->getChildren(i).get(),
                                                               minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, currentNode->getSplit(0, 0), currentPartition.getMax(), currentNode->getSplit(1, 0), currentNode->getSplit(1, 1)),
                                                               maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, currentNode->getSplit(0, 0), currentPartition.getMax(), currentNode->getSplit(1, 0), currentNode->getSplit(1, 1))));

                            }

                        }

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

                    if (result.size() % 5 == 0 && result.size() != k)
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

            memcpy(data.get() + offset, &branchFactor, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &levelsPerNode, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &fanout, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numMaxSplits, sizeof(size_t));
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

            memcpy(&branchFactor, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&levelsPerNode, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&fanout, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numMaxSplits, _data.get() + offset, sizeof(size_t));
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
            std::unique_ptr<naryTree::NodeNAry> tree = deserializeTreeRecursive(ss);

            if (root != nullptr)
                root.reset();

            buildTree(root, tree);

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            this->serializedTree = serializeTreeRecursive(root);

            ans += sizeof(size_t) + gervLib::index::Index<O, T>::getSerializedSize();
            ans += sizeof(size_t) + serializedTree.size();
            ans += sizeof(size_t) * 10;

            return ans;

        }

    };


}

#endif //GERVLIB_MVPTREE_H
