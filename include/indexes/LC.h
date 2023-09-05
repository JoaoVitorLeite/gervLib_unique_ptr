//
// Created by joaoleite on 9/5/23.
//

#ifndef GERVLIB_LC_H
#define GERVLIB_LC_H

#include "Index.h"
#include "IndexFactory.h"

namespace gervLib::index::lc
{

    template <typename O, typename T>
    class Node : serialize::Serialize {

    protected:
        std::unique_ptr<Node<O, T>> left, right;
        size_t nodeID{};

    public:
        Node() = default;

        virtual ~Node() = default;

        virtual bool isLeafNode()
        {
            return false;
        }

        void setLeft(std::unique_ptr<Node<O, T>> _left)
        {
            left.reset();
            left = std::move(_left);
        }

        void setRight(std::unique_ptr<Node<O, T>> _right)
        {
            right.reset();
            right = std::move(_right);
        }

        void setNodeID(size_t _nodeID)
        {
            nodeID = _nodeID;
        }

        std::unique_ptr<Node<O, T>>& getLeft()
        {
            return left;
        }

        std::unique_ptr<Node<O, T>>& getRight()
        {
            return right;
        }

        size_t getNodeID()
        {
            return nodeID;
        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());

            memcpy(data.get(), &nodeID, sizeof(size_t));

            return data;
        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {
            memcpy(&nodeID, _data.get(), sizeof(size_t));
            _data.reset();
        }

        size_t getSerializedSize() override
        {
            return sizeof(size_t);
        }

        virtual void clear() { }

        virtual void print(std::ostream& os) const
        {
            os << "Node type: Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;
        }

    };

    template <typename O, typename T>
    class DirectoryNode : public Node<O, T> {

    public:
        DirectoryNode() = default;

        ~DirectoryNode() override = default;

        virtual bool isLeafNode()
        {
            return false;
        }

        virtual void print(std::ostream& os) const
        {
            os << "Node type: Directory Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;
        }

    };

    template <typename O, typename T>
    class LeafNode : public Node<O, T> {

    private:
        std::unique_ptr<dataset::Dataset<O, T>> dataset;
        std::unique_ptr<index::Index<O, T>> index;
        std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot;
        double range{};

    public:
        LeafNode() = default;

        ~LeafNode() override
        {
            if (dataset != nullptr)
            {
                dataset->clear();
                dataset.reset();
            }

            if (index != nullptr) {
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

            dataset = std::move(_dataset);
        }

        void setPivot(std::unique_ptr<dataset::BasicArrayObject<O, T>> _pivot)
        {
            if (pivot != nullptr)
            {
                pivot->clear();
                pivot.reset();
            }

            pivot = std::move(_pivot);
        }

        void setIndex(std::unique_ptr<index::Index<O, T>> _index)
        {
            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }

            index = std::move(_index);
        }

        void setRange(double _range)
        {
            range = _range;
        }

        std::unique_ptr<dataset::Dataset<O, T>>& getDataset()
        {
            return dataset;
        }

        std::unique_ptr<index::Index<O, T>>& getIndex()
        {
            return index;
        }

        double getRange()
        {
            return range;
        }

        virtual bool isLeafNode()
        {
            return true;
        }

        std::unique_ptr<dataset::BasicArrayObject<O, T>>& getPivot()
        {
            return pivot;
        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &this->nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &range, sizeof(double));
            offset += sizeof(double);

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

                std::string aux = index::indexTypeMap[index->getIndexType()];
                size_t sz2 = aux.size();

                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, &sz2, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, aux.c_str(), sz2);
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

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&this->nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&range, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {

                if (dataset != nullptr)
                {
                    dataset->clear();
                    dataset.reset();
                }

                dataset = std::make_unique<dataset::Dataset<O, T>>();
                std::unique_ptr<u_char[]> datasetData = std::make_unique<u_char[]>(sz);
                memcpy(datasetData.get(), _data.get() + offset, sz);
                offset += sz;
                dataset->deserialize(std::move(datasetData));
            }
            else
                dataset = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                size_t sz2;
                std::string aux;

                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                aux.resize(sz2);
                memcpy(&aux[0], _data.get() + offset, sz2);
                offset += sz2;

                if (index != nullptr)
                {
                    index->clear();
                    index.reset();
                }

                index = index::IndexFactory<O, T>::createIndex(index::indexTypeMapReverse[aux]);

                std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
                memcpy(indexData.get(), _data.get() + offset, sz);
                offset += sz;
                index->deserialize(std::move(indexData));
            }
            else
                index = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                if (pivot != nullptr)
                {
                    pivot->clear();
                    pivot.reset();
                }

                pivot = std::make_unique<dataset::BasicArrayObject<O, T>>();
                std::unique_ptr<u_char[]> pivotData = std::make_unique<u_char[]>(sz);
                memcpy(pivotData.get(), _data.get() + offset, sz);
                offset += sz;
                pivot->deserialize(std::move(pivotData));
            }
            else
                pivot = nullptr;

            _data.reset();
        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;
            ans += sizeof(size_t);
            ans += sizeof(double);
            ans += sizeof(size_t) + (dataset != nullptr ? dataset->getSerializedSize() : 0);

            if (index != nullptr)
                ans += sizeof(size_t) * 2 + index::indexTypeMap[index->getIndexType()].size() + index->getSerializedSize();
            else
                ans += sizeof(size_t);

            ans += sizeof(size_t) + (pivot != nullptr ? pivot->getSerializedSize() : 0);

            return ans;

        }

        virtual void clear() { }

        virtual void print(std::ostream& os) const
        {
            os << "Node type: Leaf Node" << std::endl;
            os << "Node ID: " << this->nodeID << std::endl;

            if (dataset != nullptr)
                os << "Dataset: " << std::endl << *dataset << std::endl;
            else
                os << "Dataset: nullptr" << std::endl;

            if (index != nullptr)
                os << "Index: " << std::endl << *index << std::endl;
            else
                os << "Index: nullptr" << std::endl;

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

        ~LC() override = default;

        //TODO implement delete, clear, print, isEqual, buildIndex, kNN, kNNIncremental, serialize, deserialize, getSerializedSize

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

            if (useLAESA)
            {
                globalPivots = pivots::PivotFactory<O, T>::createPivot(this->pivots->getPivotType(), this->pivots);
                globalPivots->operator()(this->dataset, this->distanceFunction, this->numPivots);
            }

            root = std::make_unique<DirectoryNode<O, T>>();

            nodeQueue.push(std::make_pair(root.get(), std::move(this->dataset)));

            while (!nodeQueue.empty())
            {

                currentNode = std::move(nodeQueue.front());
                nodeQueue.pop();

                if (currentNode.second->getCardinality() <= numPerLeaf)
                {
                    std::unique_ptr<LeafNode<O, T>> leafNode = std::make_unique<LeafNode<O, T>>();
                    leafNode->setNodeID(currentNodeID++);

                    this->pivots->operator()(currentNode.second, this->distanceFunction, 1);

                    if (!storePivotsInLeaf)
                        currentNode.second->erase(this->pivots->getPivot(0));

                    double _range = std::numeric_limits<double>::min();

                    for (size_t i = 0; i < currentNode.second->getCardinality(); i++)
                    {
                        dist = this->distanceFunction->operator()(currentNode.second->getElement(i), this->pivots->getPivot(0));
                        _range = std::max(dist, _range);
                    }

                    std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot = std::make_unique<dataset::BasicArrayObject<O, T>>(this->pivots->getPivot(0));
                    leafNode->setPivot(std::move(pivot));
                    leafNode->setRange(_range);

                    if (useLAESA)
                    {
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

                    currentNode.first->setLeft(std::move(leafNode));

                }
                else
                {

                    currentNode.first->setNodeID(currentNodeID++);

                    std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>();
                    std::unique_ptr<dataset::Dataset<O, T>> rightDataset = std::make_unique<dataset::Dataset<O, T>>();
                    leftDataset->setPath(currentNode.second->getPath());
                    rightDataset->setPath(currentNode.second->getPath());
                    leftDataset->setDimensionality(currentNode.second->getDimensionality());
                    rightDataset->setDimensionality(currentNode.second->getDimensionality());
                    leftDataset->setSeed(currentNode.second->getSeed());
                    rightDataset->setSeed(currentNode.second->getSeed());

                    this->pivots->operator()(currentNode.second, this->distanceFunction, 1);

                    if (!storePivotsInLeaf)
                        currentNode.second->erase(this->pivots->getPivot(0));

                    for (size_t i = 0; i < currentNode.second->getCardinality(); i++)
                    {
                        dist = this->distanceFunction->operator()(currentNode.second->getElement(i), this->pivots->getPivot(0));
                        result->push(query::ResultEntry<dataset::BasicArrayObject<O, T>>(currentNode.second->getElement(i), dist));
                    }

                    std::vector<size_t> leftOIDs;
                    double _range = result->top().getDistance();

                    std::unique_ptr<dataset::BasicArrayObject<O, T>> pivot = std::make_unique<dataset::BasicArrayObject<O, T>>(this->pivots->getPivot(0));

                    while (!result->empty())
                    {
                        entry = result->top();
                        result->pop();
                        leftOIDs.push_back(entry.getElement().getOID());
                        leftDataset->insert(entry.getElement());
                    }

                    for (size_t i = 0; i < currentNode.second->getCardinality(); i++)
                    {
                        if (std::find(leftOIDs.begin(), leftOIDs.end(), currentNode.second->getElement(i).getOID()) == leftOIDs.end() && currentNode.second->getElement(i).getOID() != pivot->getOID())
                            rightDataset->insert(currentNode.second->getElement(i));
                    }

                    leftOIDs.clear();
                    currentNode.second->clear();
                    currentNode.second.reset();
                    result->clear();

                    if (leftDataset->getCardinality() > 0)
                    {

                        std::unique_ptr<LeafNode<O, T>> leafNode = std::make_unique<LeafNode<O, T>>();
                        leafNode->setNodeID(currentNodeID++);
                        leafNode->setPivot(std::move(pivot));
                        leafNode->setRange(_range);

                        if (useLAESA)
                        {
                            std::filesystem::path leafIndexPath(this->indexFolder);
                            leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                    this->distanceFunction->getDistanceType());
                            size_t oldDistCount = df->getDistanceCount();
                            std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                    std::move(leftDataset), std::move(df), pivots::PivotFactory<O, T>::clone(globalPivots), this->numPivots,
                                    leafIndexPath);
                            idx->getDistanceFunction()->setDistanceCount(idx->getDistanceFunction()->getDistanceCount() + oldDistCount);
                            leafNode->setIndex(std::move(idx));
                        }
                        else
                            leafNode->setDataset(std::move(leftDataset));

                        currentNode.first->setLeft(std::move(leafNode));

                    }

                    if (rightDataset->getCardinality() > 0)
                    {

                        std::unique_ptr<DirectoryNode<O, T>> directoryNode = std::make_unique<DirectoryNode<O, T>>();
                        currentNode.first->setRight(std::move(directoryNode));

                        nodeQueue.push(std::make_pair(currentNode.first->getRight().get(), std::move(rightDataset)));

                    }

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

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults, bool saveStatistics) override
        {
            throw std::runtime_error("Not implemented yet");
        }

    };

}

#endif //GERVLIB_LC_H
