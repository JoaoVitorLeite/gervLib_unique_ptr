//
// Created by joaoleite on 9/1/23.
//

#ifndef GERVLIB_PMTREE_H
#define GERVLIB_PMTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "NAryTree.h"
#include <map>

using namespace std;

namespace gervLib::index::pmtree {

    template<typename O, typename T>
    class PM_Node : public serialize::Serialize {

    public:

        PM_Node *parent_node;
        size_t node_category{}; //0: represent routing node, 1: leaf node, 2: data entry
        dataset::BasicArrayObject<O, T> feature_val;
        double dist_to_parent{};
        //size_t level{};
        MEMORY_STATUS memoryStatus = MEMORY_STATUS::NONE;

        size_t id{};
        std::vector<double> pivot_distance;

        double range{};
        std::vector<PM_Node<O, T> *> ptr_sub_tree;
        std::vector<std::pair<double, double>> hyper_rings;

//        size_t pageID{};

        std::unique_ptr<Index<O, T>> index;

        PM_Node() = default;

        PM_Node(PM_Node *parent_node_, size_t node_category_, double dist_to_parent_, size_t id_) {

            parent_node = parent_node_;
            node_category = node_category_;
            dist_to_parent = dist_to_parent_;
            id = id_;

        }

        virtual ~PM_Node()
        {
            feature_val.clear();
            pivot_distance.clear();
            hyper_rings.clear();

            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }
        }

        void clear()
        {
            //feature_val.clear();
            pivot_distance.clear();
            hyper_rings.clear();

            if (index != nullptr)
            {
                index->clear();
                index.reset();
            }
        }

        friend std::ostream &operator<<(std::ostream &os, const PM_Node &node) {
            os << "Node id: " << node.id << std::endl;
            os << "Node category: " << node.node_category << std::endl;
            os << "Feature value: " << node.feature_val << std::endl;
            os << "Range: " << node.range << std::endl;
            os << "Distance to parent: " << node.dist_to_parent << std::endl;
            os << "Memory status: " << index::memoryStatusMap[node.memoryStatus] << std::endl;

            os << "Pivot distance: ";
            for (size_t i = 0; i < node.pivot_distance.size(); i++)
                os << node.pivot_distance[i] << " ";
            os << std::endl;

            os << "Hyper rings: ";
            for (size_t i = 0; i < node.hyper_rings.size(); i++)
                os << "(" << node.hyper_rings[i].first << ", " << node.hyper_rings[i].second << ")\n";

            if (node.index != nullptr)
                os << *node.index;
            else
                os << "Index: Null\n";

            return os;

        }

        bool isEqual(PM_Node<O, T> *node)
        {

            if (memoryStatus == MEMORY_STATUS::IN_DISK || node->memoryStatus == MEMORY_STATUS::IN_DISK)
                return id == node->id;
            else
            {
                if (!feature_val.isEqual(node->feature_val))
                    return false;

                if (node_category != node->node_category)
                    return false;

                if (dist_to_parent != node->dist_to_parent)
                    return false;

                if (range != node->range)
                    return false;

                if (pivot_distance.size() != node->pivot_distance.size())
                    return false;

                for (size_t i = 0; i < pivot_distance.size(); i++)
                {
                    if (pivot_distance[i] != node->pivot_distance[i])
                        return false;
                }

                if (hyper_rings.size() != node->hyper_rings.size())
                    return false;

                for (size_t i = 0; i < hyper_rings.size(); i++)
                {
                    if (hyper_rings[i].first != node->hyper_rings[i].first || hyper_rings[i].second != node->hyper_rings[i].second)
                        return false;
                }

                if ((index == nullptr && node->index != nullptr) || (index != nullptr && node->index == nullptr))
                    return false;
                else if (index != nullptr && node->index != nullptr)
                {
                    if (!index->isEqual(node->index))
                        return false;
                }

                return true;

            }

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &node_category, sizeof(size_t));
            offset += sizeof(size_t);

            sz = feature_val.getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> feature_val_data = feature_val.serialize();
            memcpy(data.get() + offset, feature_val_data.get(), sz);
            offset += sz;
            feature_val_data.reset();

            memcpy(data.get() + offset, &dist_to_parent, sizeof(double));
            offset += sizeof(double);

            std::string aux = index::memoryStatusMap[memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;
            aux.clear();

            memcpy(data.get() + offset, &id, sizeof(size_t));
            offset += sizeof(size_t);

            sz = pivot_distance.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, pivot_distance.data(), sz * sizeof(double));
            offset += sz * sizeof(double);

            memcpy(data.get() + offset, &range, sizeof(double));
            offset += sizeof(double);

            sz = hyper_rings.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            for (size_t i = 0; i < sz; i++)
            {
                memcpy(data.get() + offset, &hyper_rings[i].first, sizeof(double));
                offset += sizeof(double);
                memcpy(data.get() + offset, &hyper_rings[i].second, sizeof(double));
                offset += sizeof(double);
            }

            if (index != nullptr)
            {
                size_t sz2;

                sz = index->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                aux = index::indexTypeMap[index->getIndexType()];
                sz2 = aux.size();
                memcpy(data.get() + offset, &sz2, sizeof(size_t));
                offset += sizeof(size_t);
                memcpy(data.get() + offset, aux.c_str(), sz2);
                offset += sz2;
                aux.clear();

                std::unique_ptr<u_char[]> index_data = index->serialize();
                memcpy(data.get() + offset, index_data.get(), sz);
                offset += sz;
                //index_data.reset();
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

            memcpy(&node_category, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> feature_val_data = std::make_unique<u_char[]>(sz);
            memcpy(feature_val_data.get(), _data.get() + offset, sz);
            offset += sz;
            feature_val.deserialize(std::move(feature_val_data));

            memcpy(&dist_to_parent, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);
            memcpy(aux.data(), _data.get() + offset, sz);
            offset += sz;
            memoryStatus = index::memoryStatusMapReverse[aux];

            memcpy(&id, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            pivot_distance.resize(sz);
            memcpy(pivot_distance.data(), _data.get() + offset, sz * sizeof(double));
            offset += sz * sizeof(double);

            memcpy(&range, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            hyper_rings.resize(sz);

            for (size_t i = 0; i < sz; i++)
            {
                memcpy(&hyper_rings[i].first, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                memcpy(&hyper_rings[i].second, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
            }

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {

                size_t sz2;

                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                aux.resize(sz2);
                memcpy(aux.data(), _data.get() + offset, sz2);
                offset += sz2;

                if (index != nullptr)
                {
                    index->clear();
                    index.reset();
                }

                index = index::IndexFactory<O, T>::createIndex(index::indexTypeMapReverse[aux]);

                std::unique_ptr<u_char[]> index_data = std::make_unique<u_char[]>(sz);
                memcpy(index_data.get(), _data.get() + offset, sz);
                offset += sz;
                index->deserialize(std::move(index_data));

            }
            else
            {
                if (index != nullptr)
                {
                    index->clear();
                    index.reset();
                }
            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = sizeof(size_t); //node category
            ans += sizeof(size_t) + feature_val.getSerializedSize(); //feature value
            ans += sizeof(double); //dist to parent
            ans += sizeof(size_t) + index::memoryStatusMap[memoryStatus].size(); //memory status
            ans += sizeof(size_t); //id
            ans += sizeof(size_t) + pivot_distance.size() * sizeof(double); //pivot distance
            ans += sizeof(double); //range
            ans += sizeof(size_t) + hyper_rings.size() * sizeof(double) * 2; //hyper rings

            if (index != nullptr)
            {
                ans += sizeof(size_t) * 2 + index::indexTypeMap[index->getIndexType()].size() + index->getSerializedSize();
            }
            else
                ans += sizeof(size_t);

            return ans;
        }

    };

    template<typename O, typename T>
    class PMTree : public Index<O, T> {

    private:
        PM_Node<O, T> *root;
        size_t numPerLeaf{}, numPivots{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::unique_ptr<pivots::Pivot<O, T>> globalPivots;
        std::string serializedTree;

    protected:
        std::string headerBuildFile() override {
            return "time,sys_time,user_time,distCount,iowrite,ioread";
        }

        std::string headerExperimentFile() override {
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

        ~PMTree()
        {
            deleteRecursive(root);
            globalPivots->clear();
            globalPivots.reset();
        }

        PM_Node<O, T>* getRoot()
        {
            return root;
        }

        void clear() override
        {
            clear_recursive(root);
        }

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;

            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);

            if (globalPivots != nullptr)
            {
                globalPivots->clear();
                globalPivots.reset();
            }

            globalPivots = pivots::PivotFactory<O, T>::clone(this->pivots);

            for (size_t i = 0; i < this->dataset->getCardinality(); i++)
            {
                insert(this->dataset->getElement(i), i);
            }

            initDisk();

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

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("PMTree::kNN not implemented yet");
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            this->leafNodeAccess = 0;
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            std::priority_queue<query::Partition<PM_Node<O, T>*>, std::vector<query::Partition<PM_Node<O, T>*>>, std::greater<query::Partition<PM_Node<O, T>*>>> nodeQueue;
            std::priority_queue<query::ResultEntry<O>, std::vector<query::ResultEntry<O>>, std::greater<query::ResultEntry<O>>> elementQueue;
            query::Result<O> result;
            result.setMaxSize(k);
            query::Partition<PM_Node<O, T>*> currentPartition;
            PM_Node<O, T>* currentNode;
            LAESA<O, T>* laesa;
            std::vector<double> query_to_pivot;
            update_pivot_distance(query_to_pivot, query);

            nodeQueue.push(query::Partition<PM_Node<O, T>*>(root, 0.0, std::numeric_limits<double>::max()));

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty())) {

                if (elementQueue.empty()) {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->memoryStatus == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->id);
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if (is_leaf_node(currentNode)) {

                        this->leafNodeAccess++;

                        if (currentNode->index == nullptr) {

                            for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++) {
                                elementQueue.push(
                                        query::ResultEntry<O>(currentNode->ptr_sub_tree[i]->feature_val.getOID(),
                                                              this->distanceFunction->operator()(
                                                                      currentNode->ptr_sub_tree[i]->feature_val,
                                                                      query)));
                            }

                        }
                        else
                        {
                            laesa = (LAESA<O, T>*) currentNode->index.get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentNode->index->getPrunning();

                        }

                    } else {
                        for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++) {

                            nodeQueue.push(query::Partition<PM_Node<O, T>*>(currentNode->ptr_sub_tree[i],
                                                              minDistNode(currentNode->ptr_sub_tree[i], query,
                                                                          query_to_pivot),
                                                              maxDistNode(currentNode->ptr_sub_tree[i], query,
                                                                          query_to_pivot)));

                        }
                    }

                    if (currentNode->memoryStatus == MEMORY_STATUS::IN_DISK)
                    {
                        currentNode->clear();
                    }

                } else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance()) {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->memoryStatus == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->id);
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if (is_leaf_node(currentNode)) {

                        this->leafNodeAccess++;

                        if (currentNode->index == nullptr) {

                            for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++) {
                                elementQueue.push(
                                        query::ResultEntry<O>(currentNode->ptr_sub_tree[i]->feature_val.getOID(),
                                                              this->distanceFunction->operator()(
                                                                      currentNode->ptr_sub_tree[i]->feature_val,
                                                                      query)));
                            }

                        }
                        else
                        {
                            laesa = (LAESA<O, T>*) currentNode->index.get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentNode->index->getPrunning();

                        }

                    } else {
                        for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++) {

                            nodeQueue.push(query::Partition<PM_Node<O, T>*>(currentNode->ptr_sub_tree[i],
                                                              minDistNode(currentNode->ptr_sub_tree[i], query,
                                                                          query_to_pivot),
                                                              maxDistNode(currentNode->ptr_sub_tree[i], query,
                                                                          query_to_pivot)));

                        }
                    }

                    if (currentNode->memoryStatus == MEMORY_STATUS::IN_DISK)
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

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<PMTree<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

        }

        void print(std::ostream& os) const override {

            std::stack<PM_Node<O, T>*> nodeStack;
            nodeStack.push(root);

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "PMTree" << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Number of objects per leaf: " << numPerLeaf << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;
            os << "Use LAESA: " << (useLAESA ? "true" : "false") << std::endl;

            while (!nodeStack.empty())
            {

                PM_Node<O, T>* currentNode = nodeStack.top();
                nodeStack.pop();

                os << "**********************************************************************************************************************************************************************************\n\n";
                os << *currentNode << std::endl;

                if (currentNode->node_category != 1)
                {
                    for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++)
                        nodeStack.push(currentNode->ptr_sub_tree[i]);
                }

            }

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &numPerLeaf, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numPivots, sizeof(size_t));
            offset += sizeof(size_t);

            sz = storeLeafNode ? 1 : 0;
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = storeDirectoryNode ? 1 : 0;
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = useLAESA ? 1 : 0;
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = gervLib::index::Index<O, T>::getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> index_data = gervLib::index::Index<O, T>::serialize();
            memcpy(data.get() + offset, index_data.get(), sz);
            offset += sz;
            index_data.reset();

            if (globalPivots != nullptr)
            {
                sz = globalPivots->getSerializedSize();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::string aux = pivots::PIVOT_TYPE2STR[globalPivots->getPivotType()];
                size_t sz2 = aux.size();

                memcpy(data.get() + offset, &sz2, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, aux.c_str(), sz2);
                offset += sz2;
                aux.clear();

                std::unique_ptr<u_char[]> pivots_data = globalPivots->serialize();
                memcpy(data.get() + offset, pivots_data.get(), sz);
                offset += sz;
                pivots_data.reset();
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            sz = serializedTree.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);
            serializedTree.clear();

            memcpy(data.get() + offset, serializedTree.c_str(), sz);
            offset += sz;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&numPerLeaf, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numPivots, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeLeafNode = sz == 1;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeDirectoryNode = sz == 1;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            useLAESA = sz == 1;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> index_data = std::make_unique<u_char[]>(sz);
            memcpy(index_data.get(), _data.get() + offset, sz);
            offset += sz;
            gervLib::index::Index<O, T>::deserialize(std::move(index_data));

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                size_t sz2;
                std::string aux;

                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                aux.resize(sz2);
                memcpy(aux.data(), _data.get() + offset, sz2);
                offset += sz2;

                if (globalPivots != nullptr)
                {
                    globalPivots->clear();
                    globalPivots.reset();
                }

                globalPivots = pivots::PivotFactory<O, T>::createPivot(pivots::STR2PIVOT_TYPE[aux]);

                std::unique_ptr<u_char[]> pivots_data = std::make_unique<u_char[]>(sz);
                memcpy(pivots_data.get(), _data.get() + offset, sz);
                offset += sz;
                globalPivots->deserialize(std::move(pivots_data));
            }
            else
            {
                globalPivots = nullptr;
            }

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            serializedTree.resize(sz);
            memcpy(serializedTree.data(), _data.get() + offset, sz);
            offset += sz;

            _data.reset();

            std::stringstream ss(serializedTree);
            std::unique_ptr<naryTree::NodeNAry> aux = deserializeTreeRecursive(ss);

            if (root != nullptr)
            {
                root->clear();
                delete root;
                root = nullptr;
            }

            //root = new PM_Node<O, T>();

            buildTree(root, aux);

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;
            this->serializedTree = serializeTreeRecursive(root);

            ans += sizeof(size_t) + gervLib::index::Index<O, T>::getSerializedSize();
            ans += sizeof(size_t) * 5;
            ans += sizeof(size_t) + (globalPivots == nullptr ? 0 : globalPivots->getSerializedSize() + sizeof(size_t) + pivots::PIVOT_TYPE2STR[globalPivots->getPivotType()].size());
            ans += sizeof(size_t) + serializedTree.size();

            return ans;
        }


    private:
        void deleteRecursive(PM_Node<O, T>* node)
        {
            if (node == nullptr)
                return;

            for (size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                deleteRecursive(node->ptr_sub_tree[i]);

            node->clear();
            delete node;
            node = nullptr;
        }

        bool isEqualHelper(PM_Node<O, T>* node1, PM_Node<O, T>* node2)
        {
            if (node1 == nullptr && node2 == nullptr)
                return true;

            if ((node1 != nullptr && node2 == nullptr) || (node1 == nullptr && node2 != nullptr))
                return false;

            if (!node1->isEqual(node2))
                return false;

            if (node1->ptr_sub_tree.size() != node2->ptr_sub_tree.size())
                return false;

            for (size_t i = 0; i < node1->ptr_sub_tree.size(); i++)
            {
                if (!isEqualHelper(node1->ptr_sub_tree[i], node2->ptr_sub_tree[i]))
                    return false;
            }

            return true;

        }

        void clear_recursive(PM_Node<O, T>* node)
        {

            if (node == nullptr)
                return;

            if (is_leaf_node(node))
            {
                node->clear();
            }
            else
            {
                for (size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                    clear_recursive(node->ptr_sub_tree[i]);
            }

        }

        void buildTree(PM_Node<O, T>*& node, std::unique_ptr<naryTree::NodeNAry>& aux)
        {

            if (aux == nullptr)
                return;

            if (aux->value[0] == 'L')
            {
                node = new PM_Node<O, T>();
                node->id = std::stoull(aux->value.substr(1));

                std::unique_ptr<u_char[]> data = this->pageManager->load(node->id);
                node->deserialize(std::move(data));

                if (!storeLeafNode) {
                    node->memoryStatus = index::MEMORY_STATUS::IN_MEMORY;
                }
                else {
                    node->memoryStatus = index::MEMORY_STATUS::IN_DISK;
                    node->clear();
                }
            }
            else
            {
                node = new PM_Node<O, T>();
                node->id = std::stoull(aux->value.substr(1));

                std::unique_ptr<u_char[]> data = this->pageManager->load(node->id);
                node->deserialize(std::move(data));

                if (!storeDirectoryNode) {
                    node->memoryStatus = index::MEMORY_STATUS::IN_MEMORY;
                }
                else {
                    node->memoryStatus = index::MEMORY_STATUS::IN_DISK;
                    node->clear();
                }
            }

            for (size_t i = 0; i < aux->children.size(); i++)
            {
                node->ptr_sub_tree.push_back(new PM_Node<O, T>());
                buildTree(node->ptr_sub_tree[i], aux->children[i]);
            }

        }

        std::string serializeTreeRecursive(PM_Node<O, T>* node)
        {
            if (node == nullptr)
                return "null ";

            if (node->memoryStatus != index::MEMORY_STATUS::IN_DISK)
            {
                std::unique_ptr<u_char[]> data = node->serialize();
                this->pageManager->save(node->id, std::move(data), node->getSerializedSize());
                node->memoryStatus = index::MEMORY_STATUS::IN_DISK;
            }

            std::string result = (is_leaf_node(node) ? "L" : "D") + std::to_string(node->id) + " " + std::to_string(node->ptr_sub_tree.size()) + " ";

            for (size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                result += serializeTreeRecursive(node->ptr_sub_tree[i]);

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

        void initDisk()
        {
            std::queue<PM_Node<O, T>*> nodeQueue;
            PM_Node<O, T>* currentNode;
            size_t currentNodeId = 0;

            nodeQueue.push(root);

            while (!nodeQueue.empty())
            {

                currentNode = nodeQueue.front();
                nodeQueue.pop();
                currentNode->id = currentNodeId++;

                if (is_leaf_node(currentNode))
                {
                    std::unique_ptr<dataset::Dataset<O, T>> laesaDataset = std::make_unique<dataset::Dataset<O, T>>();
                    std::unique_ptr<matrix::Matrix<O, T>> matrix = std::make_unique<matrix::Matrix<O, T>>(this->numPivots, currentNode->ptr_sub_tree.size());

                    for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++)
                    {
                        laesaDataset->insert(currentNode->ptr_sub_tree[i]->feature_val);
                        currentNode->ptr_sub_tree[i]->feature_val.clear();

                        for (size_t j = 0; j < this->numPivots; j++)
                            matrix->setValue(j, i, currentNode->ptr_sub_tree[i]->pivot_distance[j]);

                        currentNode->ptr_sub_tree[i]->clear();
                        delete currentNode->ptr_sub_tree[i];
                        currentNode->ptr_sub_tree[i] = nullptr;

                    }

                    currentNode->ptr_sub_tree.clear();

                    std::filesystem::path leafIndexPath(this->indexFolder);
                    leafIndexPath /= "laesa_leafnode_" + std::to_string(currentNode->id);
                    std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                            this->distanceFunction->getDistanceType());
                    std::unique_ptr<LAESA<O, T>> idx = std::make_unique<index::LAESA<O, T>>();
                    idx->setDataset(std::move(laesaDataset));
                    idx->setDistanceFunction(std::move(df));
                    idx->setPivots(pivots::PivotFactory<O, T>::clone(globalPivots));
                    idx->setIndexFolder(leafIndexPath);
                    idx->setMatrix(std::move(matrix));

                    if (currentNode->index != nullptr) {
                        currentNode->index->clear();
                        currentNode->index.reset();
                    }

                    currentNode->index = std::move(idx);

                    if (storeLeafNode)
                    {
                        currentNode->memoryStatus = gervLib::index::MEMORY_STATUS::IN_DISK;
                        std::unique_ptr<u_char[]> leafData = currentNode->serialize();
                        this->pageManager->save(currentNode->id, std::move(leafData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->memoryStatus = gervLib::index::MEMORY_STATUS::IN_MEMORY;

                }
                else
                {

                    for (size_t i = 0; i < currentNode->ptr_sub_tree.size(); i++)
                        nodeQueue.push(currentNode->ptr_sub_tree[i]);

                    if (storeDirectoryNode) {
                        currentNode->memoryStatus = gervLib::index::MEMORY_STATUS::IN_DISK;
                        std::unique_ptr<u_char[]> directoryData = currentNode->serialize();
                        this->pageManager->save(currentNode->id, std::move(directoryData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->memoryStatus = gervLib::index::MEMORY_STATUS::IN_MEMORY;

                }

            }
        }

        void insert(dataset::BasicArrayObject<O, T>& feature_val_, size_t id_)
        {
            auto* new_node = new PM_Node<O, T>(nullptr, 2, -1, id_);
            new_node->feature_val = feature_val_;
            new_node->range = -1;
            new_node->id = id_; //DUVIDA
            update_pivot_distance(new_node->pivot_distance, feature_val_);

            if(root == nullptr)
            {

                auto* new_root = new PM_Node<O, T>(nullptr, 1, -1, -1);
                new_root->feature_val = feature_val_;
                new_root->ptr_sub_tree.push_back(new_node);
                new_root->range = cal_cover_radius(new_root);
                update_hyper_rings(new_root->hyper_rings, feature_val_);
                root = new_root;

                new_node->parent_node = new_root;
                new_node->dist_to_parent = cal_dist_to_parent(new_node);

            }
            else
            {

                insert(&root, &new_node);

            }

        }

        void insert(PM_Node<O, T> ** cur_node_address_, PM_Node<O, T> ** insert_node_address_)
        {

            PM_Node<O, T>* cur_node_ = *cur_node_address_;
            PM_Node<O, T>* insert_node_ = *insert_node_address_;

            if(is_leaf_node(cur_node_))
            {

                if(is_full(cur_node_))
                {

                    split(&cur_node_, &insert_node_);

                }
                else
                {

                    cur_node_->range = cal_cover_radius(cur_node_, insert_node_);
                    cur_node_->ptr_sub_tree.push_back(insert_node_);

                    insert_node_->parent_node = cur_node_;
                    insert_node_->dist_to_parent = cal_dist_to_parent(insert_node_);
                    update_hyper_rings(cur_node_->hyper_rings, insert_node_->feature_val);

                    //MUDANCAS JOAO
                    PM_Node<O, T>* update_ptr = cur_node_->parent_node;

                    while(update_ptr != nullptr)
                    {

                        update_ptr->range = cal_cover_radius(update_ptr);
                        merge_subNode_HR(update_ptr);
                        update_ptr = update_ptr->parent_node;

                    }
                    //MUDANCAS JOAO

                }

            }
            else
            {

                PM_Node<O, T>** next_node = get_next_node_and_update_range(&cur_node_, &insert_node_);
                update_hyper_rings((*next_node)->hyper_rings, insert_node_->feature_val);
                insert(next_node, &insert_node_);

            }

        }

        PM_Node<O, T>** get_next_node_and_update_range(PM_Node<O, T> ** cur_node_address_, PM_Node<O, T> ** insert_node_address_)
        {

            PM_Node<O, T>* cur_node_ = *cur_node_address_;
            PM_Node<O, T>* insert_node_ = *insert_node_address_;

            std::vector<std::pair<double, PM_Node<O, T>**>> data_vec;
            std::vector<std::pair<double, PM_Node<O, T>**>> data2_vec;

            for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
            {

                double dist = this->distanceFunction->operator()(cur_node_->ptr_sub_tree[i]->feature_val, insert_node_->feature_val);
                //std::cout << "range = " << cur_node_->ptr_sub_tree[i]->range << "\n";

                if(dist <= cur_node_->ptr_sub_tree[i]->range)
                {

                    //std::cout << "DIS in range= "  << cur_node_->ptr_sub_tree[i]->range - dist << "\n";
                    data_vec.push_back(std::make_pair(dist, &(cur_node_->ptr_sub_tree[i])));

                }

                if(data_vec.size() == 0)
                {

                    dist -= cur_node_->ptr_sub_tree[i]->range;
                    //std::cout << "DIS out range= "  << dist << "\n";
                    data2_vec.push_back(std::make_pair(dist, &(cur_node_->ptr_sub_tree[i])));

                }

            }

            if(data_vec.size() != 0)
            {

                std::sort(data_vec.begin(), data_vec.end());
                return data_vec[0].second;

            }

            std::sort(data2_vec.begin(), data2_vec.end());
            (*(data2_vec[0].second))->range = data2_vec[0].first + (*(data2_vec[0].second))->range;
            //std::cout << (*data2_vec[0].second)->feature_val.getOID() << "\n";
            return data2_vec[0].second;

        }

        void update_pivot_distance(std::vector<double>& pivot_distance_, dataset::BasicArrayObject<O, T>& data_feature_val_)
        {

            if(pivot_distance_.empty())
            {

                pivot_distance_.resize(this->numPivots, 0);

            }

            for(size_t i = 0; i < this->numPivots; ++i)
            {

                double dist = this->distanceFunction->operator()(data_feature_val_, globalPivots->getPivot(i));
                pivot_distance_[i] = dist;

            }

        }

        double cal_cover_radius(const PM_Node<O, T>* node_)
        {

            if(node_->ptr_sub_tree.size() == 0)
            {

                return 0;

            }

            double max_range = 0.0;

            if(is_leaf_node(node_))
            {

                for(size_t i = 0; i < node_->ptr_sub_tree.size(); ++i)
                {

                    double dist = this->distanceFunction->operator()(node_->ptr_sub_tree[i]->feature_val, node_->feature_val);

                    if(dist > max_range)
                    {

                        max_range = dist;

                    }

                }

            }
            else
            {

                for(size_t i = 0; i < node_->ptr_sub_tree.size(); ++i)
                {

                    double dist = this->distanceFunction->operator()(node_->ptr_sub_tree[i]->feature_val, node_->feature_val);
                    dist += node_->ptr_sub_tree[i]->range;

                    if(dist > max_range)
                    {

                        max_range = dist;

                    }

                }

            }

            return max_range;

        }

        double cal_dist_to_parent(const PM_Node<O, T>* node_)
        {

            double dist;

            if(!(node_ == nullptr || node_->parent_node == nullptr))
            {

                dist = this->distanceFunction->operator()(node_->feature_val, node_->parent_node->feature_val);

            }
            else
            {

                throw std::invalid_argument("Error in cal_dist_to_parent \n");

            }

            return dist;

        }

        void update_hyper_rings(std::vector<std::pair<double,double>>& hyper_rings_, dataset::BasicArrayObject<O, T>& data_feature_val_)
        {

            if(hyper_rings_.empty())
            {

                hyper_rings_.resize(this->numPivots, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::min()));

            }

            for(size_t i = 0; i < this->numPivots; ++i)
            {

                double dist = this->distanceFunction->operator()(data_feature_val_, globalPivots->getPivot(i));

                if(dist < hyper_rings_[i].first)
                {

                    hyper_rings_[i].first = dist;

                }

                if(dist > hyper_rings_[i].second)
                {

                    hyper_rings_[i].second = dist;

                }
            }

        }

        bool is_leaf_node(const PM_Node<O, T>* node_)
        {

            return node_->node_category == 1;

        }

        bool is_full(PM_Node<O, T>* node_)
        {

            return node_->ptr_sub_tree.size() >= this->numPerLeaf;

        }

        double cal_cover_radius(const PM_Node<O, T>* node_, const PM_Node<O, T>* cmp_node_)
        {

            double dist;

            if(is_leaf_node(node_))
            {

                dist = this->distanceFunction->operator()(cmp_node_->feature_val, node_->feature_val);

            }
            else
            {

                dist = this->distanceFunction->operator()(cmp_node_->feature_val, node_->feature_val);
                dist += cmp_node_->range;

            }

            if(node_->range < dist)
            {

                return dist;

            }
            else
            {

                return node_->range;

            }

        }

        void merge_subNode_HR(PM_Node<O, T>* cur_node_)
        {

            if(is_leaf_node(cur_node_))
            {

                for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
                {

                    update_hyper_rings(cur_node_->hyper_rings, cur_node_->ptr_sub_tree[i]->feature_val);

                }

            }
            else
            {

                if(cur_node_->hyper_rings.size() == 0)
                {

                    cur_node_->hyper_rings.resize(this->numPivots, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::min()));

                }

                for(size_t i = 0; i < this->numPivots; i++)
                {

                    for(size_t j = 0; j < cur_node_->ptr_sub_tree.size(); ++j)
                    {

                        double sub_min;
                        double sub_max;

                        if(cur_node_->ptr_sub_tree[j]->ptr_sub_tree.size() != 0)
                        {

                            sub_min = cur_node_->ptr_sub_tree[j]->hyper_rings[i].first;
                            sub_max = cur_node_->ptr_sub_tree[j]->hyper_rings[i].second;

                        }
                        else
                            continue;
                        double cur_min = cur_node_->hyper_rings[i].first;
                        double cur_max = cur_node_->hyper_rings[i].second;

                        if(cur_min > sub_min)
                        {

                            cur_node_->hyper_rings[i].first = cur_node_->ptr_sub_tree[j]->hyper_rings[i].first;

                        }

                        if(cur_max < sub_max)
                        {

                            cur_node_->hyper_rings[i].second = cur_node_->ptr_sub_tree[j]->hyper_rings[i].second;

                        }

                    }

                }

            }

        }

        void split(PM_Node<O, T> ** cur_node_ptr_address_, PM_Node<O, T> ** insert_node_address_)
        {

            //std::cout << "SPLIT\n";
            PM_Node<O, T> * cur_node_ = *cur_node_ptr_address_;
            PM_Node<O, T> * insert_node_ = *insert_node_address_;

            PM_Node<O, T>* new_Mnode_1 = new PM_Node<O, T>(nullptr, 0, -1, -1);
            PM_Node<O, T>* new_Mnode_2 = new PM_Node<O, T>(nullptr, 0, -1, -1);

            std::vector<PM_Node<O, T>*> entries = cur_node_->ptr_sub_tree;
            entries.push_back(insert_node_);

            promote(entries, new_Mnode_1, new_Mnode_2);

            partition(entries, new_Mnode_1, new_Mnode_2);

            if(is_root_node(cur_node_))
            {

                PM_Node<O, T>* new_root = new PM_Node<O, T>(nullptr, 0, -1, -1);
                new_root->feature_val = new_Mnode_1->feature_val;

                new_Mnode_1->parent_node = new_root;
                new_Mnode_1->dist_to_parent = cal_dist_to_parent(new_Mnode_1);
                new_Mnode_2->parent_node = new_root;
                new_Mnode_2->dist_to_parent = cal_dist_to_parent(new_Mnode_2);
                merge_subNode_HR(new_Mnode_1);
                merge_subNode_HR(new_Mnode_2);

                new_root->ptr_sub_tree.push_back(new_Mnode_1);
                new_root->ptr_sub_tree.push_back(new_Mnode_2);
                new_root->range = cal_cover_radius(new_root);
                merge_subNode_HR(new_root);
                root = new_root;

            }
            else
            {

                new_Mnode_1->parent_node = cur_node_->parent_node;
                new_Mnode_1->dist_to_parent = cal_dist_to_parent(new_Mnode_1);
                new_Mnode_2->parent_node = cur_node_->parent_node;
                new_Mnode_2->dist_to_parent = cal_dist_to_parent(new_Mnode_2);
                merge_subNode_HR(new_Mnode_1); //JOAO
                assign_node_all_value(cur_node_, new_Mnode_1);
                //merge_subNode_HR(cur_node_); //JOAO: PROVAVEL COMANDO SEM UTILIDADE
                merge_subNode_HR(new_Mnode_2);

                delete (new_Mnode_1);

                if(is_full(new_Mnode_2->parent_node))
                {

                    split(&(new_Mnode_2->parent_node), &new_Mnode_2);

                }
                else
                {

                    new_Mnode_2->parent_node->ptr_sub_tree.push_back(new_Mnode_2);

                }

                //MUDANCAS JOAO
                PM_Node<O, T>* update_ptr = cur_node_->parent_node;

                while(update_ptr != nullptr)
                {

                    update_ptr->range = cal_cover_radius(update_ptr);
                    merge_subNode_HR(update_ptr);
                    update_ptr = update_ptr->parent_node;

                }
                //MUDANCAS JOAO


//        std::cout << "DPS 1 = " << root->range << std::endl;
//        root->range = cal_cover_radius(root);
//        std::cout << "DPS 2 = " << root->range << "\n\n";

            }

        }

        void assign_node_all_value(PM_Node<O, T>* cur_node_, PM_Node<O, T>* new_node_)
        {

            cur_node_->dist_to_parent = new_node_->dist_to_parent;
            cur_node_->feature_val = new_node_->feature_val;
            cur_node_->id = new_node_->id;
            cur_node_->node_category = new_node_->node_category;
            cur_node_->parent_node = new_node_->parent_node;
            cur_node_->ptr_sub_tree = new_node_->ptr_sub_tree;
            cur_node_->range = new_node_->range;
            cur_node_->pivot_distance = new_node_->pivot_distance;
            //cur_node_->hyper_rings = cur_node_->hyper_rings; //DUVIDA
            cur_node_->hyper_rings = new_node_->hyper_rings; //DUVIDA

            for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
            {

                cur_node_->ptr_sub_tree[i]->parent_node = cur_node_;

            }

        }

        bool is_root_node(const PM_Node<O, T>* node_)
        {

            return node_->parent_node == nullptr;

        }

        void promote(std::vector<PM_Node<O, T>*>& entries_, PM_Node<O, T>* node_1_, PM_Node<O, T>* node_2_)
        {

            std::unique_ptr<dataset::Dataset<O, T>> data = std::make_unique<dataset::Dataset<O, T>>();

            for (size_t i = 0; i < entries_.size(); i++)
            {
                data->insert(entries_[i]->feature_val);
            }

//            std::vector<dataset::BasicArrayObject<O, T>> v;
//
//            for(size_t i = 0; i < entries_.size(); i++)
//            {
//
//                v.push_back(entries_[i]->feature_val);
//
//            }
//
//            dataset::Dataset<O, T>* data = new dataset::Dataset<O, T>(v, entries_.size(), entries_[0]->feature_val.size());

            this->pivots->operator()(data, this->distanceFunction, 2);
//            pvt->generatePivots(data, df, 2);

            node_1_->feature_val = this->pivots->getPivot(0);
            node_2_->feature_val = this->pivots->getPivot(1);

            //std::cout << "PROMOTE : " << node_1_->feature_val.toStringWithOID() << "\n";
            //std::cout << "PROMOTE : " << node_2_->feature_val.toStringWithOID() << "\n\n";

//            v.clear();
//            delete (data);
            data->clear();
            data.reset();

        }

        void partition(std::vector<PM_Node<O, T>*>& entries_, PM_Node<O, T>* node_1_, PM_Node<O, T>* node_2_)
        {

            if(entries_[0]->node_category == 1)
            {

                node_1_->node_category = node_2_->node_category = 0;

            }
            else if(entries_[0]->node_category == 2)
            {

                node_1_->node_category = node_2_->node_category = 1;

            }
            else if(entries_[0]->node_category == 0)
            {

                node_1_->node_category = node_2_->node_category = 0;

            }
            else
            {

                throw std::invalid_argument("Error in partition");

            }

            for(size_t i = 0; i < entries_.size(); i++)
            {

                if(this->distanceFunction->operator()(entries_[i]->feature_val, node_1_->feature_val) <= this->distanceFunction->operator()(entries_[i]->feature_val, node_2_->feature_val))
                {

                    node_1_->ptr_sub_tree.push_back(entries_[i]);
                    entries_[i]->parent_node = node_1_;
                    //std::cout << "NODE 1 - PARTITION\n";

                }
                else
                {

                    node_2_->ptr_sub_tree.push_back(entries_[i]);
                    entries_[i]->parent_node = node_2_;
                    //std::cout << "NODE 2 - PARTITION\n";

                }

                entries_[i]->dist_to_parent = cal_dist_to_parent(entries_[i]);

            }

            node_1_->range = cal_cover_radius(node_1_);
            node_2_->range = cal_cover_radius(node_2_);

        }

        double minDistNode(PM_Node<O, T>* cur_node_, dataset::BasicArrayObject<O, T>& element, std::vector<double> dist_to_query)
        {

//    double dist = df->getDistance(cur_node_->feature_val, element), ans;
//    std::cout << cur_node_->feature_val.toStringWithOID() << "\n";
//    std::cout << element.toStringWithOID() << "\n";
//    std::cout << "MINDIST = " << dist << " / RANGE = " << cur_node_->range << "\n\n\n";

//    if(dist <= cur_node_->range) //Dentro da bola
//    {

//        ans = 0.0;

//    }
//    else //Fora da bola
//    {

//        ans = dist - cur_node_->range;

//    }

//    return ans;

            double d_hr_max = 0.0, d_hr_min = 0.0, dist_p_i_q;

            for(size_t i = 0; i < this->numPivots; i++)
            {

//        dist_p_i_q = df->getDistance(pivot_vec.getFeatureVector(i), element);
                dist_p_i_q = dist_to_query[i];
                d_hr_max = std::max(d_hr_max, dist_p_i_q - cur_node_->hyper_rings[i].second);
                d_hr_min = std::max(d_hr_min, cur_node_->hyper_rings[i].first - dist_p_i_q);

            }

            return std::max({0.0, this->distanceFunction->operator()(cur_node_->feature_val, element) - cur_node_->range, d_hr_max, d_hr_min});

        }

        double maxDistNode(PM_Node<O, T>* cur_node_, dataset::BasicArrayObject<O, T>& element, std::vector<double> dist_to_query)
        {

            //return df->getDistance(cur_node_->feature_val, element) + cur_node_->range;

            double d_hr = std::numeric_limits<double>::max();

            for(size_t i = 0; i < this->numPivots; i++)
            {

//        d_hr = std::min(d_hr, df->getDistance(pivot_vec.getFeatureVector(i), element) + cur_node_->hyper_rings[i].second);
                d_hr = std::min(d_hr, dist_to_query[i] + cur_node_->hyper_rings[i].second);

            }

            return std::min(this->distanceFunction->operator()(cur_node_->feature_val, element) + cur_node_->range, d_hr);

        }

    };

}

#endif //GERVLIB_PMTREE_H
