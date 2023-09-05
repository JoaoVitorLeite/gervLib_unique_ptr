//
// Created by joaoleite on 9/4/23.
//

#ifndef GERVLIB_SPBTREE_H
#define GERVLIB_SPBTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "HilbertCurve.h"
#include "EquiDepth.h"
#include "btree_multimap.h"

namespace gervLib::index::spbtree
{

    typedef stx::btree_multimap<size_t, double, unsigned long long, std::unique_ptr<dataset::BasicArrayObject<size_t, double>>, std::less<> > btree_type;

    template <typename O, typename T>
    class SPBTree : public Index<O, T>
    {

    private:
        std::unique_ptr<hilbert::HilbertCurve> hc;
        std::unique_ptr<equidepth::EquiDepth<double>> ed;
        size_t numPerLeaf{}, numPivots{}, numBins{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::unique_ptr<stx::btree_multimap<O, T, unsigned long long, std::unique_ptr<dataset::BasicArrayObject<size_t, double>>, std::less<> >> bptree;

    private:
        std::unique_ptr<std::vector<unsigned long long>> keys_serialize;

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
        SPBTree()
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->bptree = nullptr;
            this->hc = nullptr;
            this->ed = nullptr;
            this->numPivots = 0;
            this->numPerLeaf = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->indexFolder = "";
        }

        SPBTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset, std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
                std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _numBins, size_t _pageSize = 0,
                bool _storeDirectoryNode = false, bool _storeLeafNode = false, bool _useLAESA = true, std::string folder="")
        {
            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->numPivots = std::min((size_t)7, _numPivots);
            this->pageSize = _pageSize;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numBins = _numBins;

            unsigned long long p = (unsigned long long)log2(numBins-1) + 1;

            this->bptree = std::make_unique<btree_type>();
            this->hc = std::make_unique<hilbert::HilbertCurve>(p, numPivots);
            this->ed = std::make_unique<equidepth::EquiDepth<T>>(numBins, numPivots);
            this->numPerLeaf = _numPerLeaf;
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("spb_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit SPBTree(std::string _folder, std::string serializedFile = "")
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->bptree = nullptr;
            this->hc = nullptr;
            this->ed = nullptr;
            this->numPivots = 0;
            this->numPerLeaf = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~SPBTree() override
        {
            ed.reset();
            hc.reset();
            this->bptree->clear();
        }

        //TODO implement isEqual, serialize, deserialize, getSerializedSize

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;

            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);

            std::vector<std::vector<double>> pivot_mapping;
            pivot_mapping.resize(this->dataset->getCardinality(), std::vector<double>(numPivots));

            std::vector<std::vector<unsigned long long>> disc;
            disc.resize(this->dataset->getCardinality(), std::vector<unsigned long long>(numPivots));

            std::vector<unsigned long long> keys(this->dataset->getCardinality());

            std::ofstream file_pivot_mapping(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_pivot_mapping.txt");
            std::ofstream file_disc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_disc.txt");
            std::ofstream file_sfc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_sfc.txt");

            for(size_t i = 0; i < this->dataset->getCardinality(); i++)
            {

                for(size_t j = 0; j < numPivots; j++)
                {

                    pivot_mapping[i][j] = this->distanceFunction->operator()(this->dataset->getElement(i), this->pivots->getPivot(j));
                    file_pivot_mapping << pivot_mapping[i][j] << " ";

                }

                file_pivot_mapping << "\n";

            }

            ed->build(pivot_mapping);
            ed->saveToFile(this->indexFolder);
            file_pivot_mapping.close();

            for(size_t i = 0; i < this->dataset->getCardinality(); i++)
            {

                for(size_t j = 0; j < numPivots; j++)
                {

                    disc[i][j] = ed->getBin(j, pivot_mapping[i][j]);
                    file_disc << disc[i][j] << " ";

                }

                file_disc << "\n";

                keys[i] = hc->distance_from_point(disc[i]);
                file_sfc << keys[i] << "\n";

            }

            pivot_mapping.clear();
            disc.clear();
            ed->clear();
            file_disc.close();
            file_sfc.close();

            std::vector<std::pair<unsigned long long, std::unique_ptr<dataset::BasicArrayObject<size_t, double>>>> insertValues(keys.size());

            for(size_t i = 0; i < keys.size(); i++)
            {

                std::unique_ptr<dataset::BasicArrayObject<size_t, double>> obj = std::make_unique<dataset::BasicArrayObject<size_t, double>>(this->dataset->getElement(i));
                insertValues[i] = std::make_pair(keys[i], std::move(obj));

            }

            keys.clear();

            std::sort(insertValues.begin(), insertValues.end());
            bptree->bulk_load(insertValues.begin(), insertValues.end());

            for (auto& insertValue : insertValues)
            {

                insertValue.second.reset();

            }

            insertValues.clear();
            this->dataset->clear();
            this->dataset.reset();

            bptree->setHilbertCurve(hc.get());
            bptree->setNumberOfPivots(numPivots);
            bptree->init_Key_Minmax();
            bptree->setPageManager(this->pageManager.get());
            bptree->setPivots(this->pivots);
            bptree->setIndexFolder(this->indexFolder);
            bptree->setDistanceFunction(this->distanceFunction.get());
            bptree->initDisk(storeDirectoryNode, storeLeafNode, useLAESA);

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

            utils::Timer timer{};
            ed->load(this->indexFolder);
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            this->leafNodeAccess = 0;
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            std::priority_queue<query::Partition<btree_type::Node*>, std::vector<query::Partition<btree_type::Node*>>, std::greater<>> nodeQueue;
            std::priority_queue<query::ResultEntry<O>, std::vector<query::ResultEntry<O>>, std::greater<query::ResultEntry<O>>> elementQueue;
            query::Result<O> result;
            result.setMaxSize(k);
            query::Partition<btree_type::Node*> currentPartition;
            btree_type::Node *currentNode;
            btree_type::Leaf *leafnode;
            btree_type::Inner *innernode;
            LAESA<O, T>* laesa;

            auto sq_ = std::vector<double>(numPivots);

            for(size_t i = 0; i < numPivots; i++)
            {

                sq_[i] = this->distanceFunction->operator()(query, this->pivots->getPivot(i));

            }

            nodeQueue.push(query::Partition<btree_type::Node*>(bptree->getRoot(), 0.0, std::numeric_limits<double>::max()));

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty())) {

                if (elementQueue.empty()) {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->memory_status == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->nodeID);
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if(currentNode->isleafnode()) {

                        leafnode = dynamic_cast<btree_type::Leaf *>(currentNode);
                        this->leafNodeAccess++;

                        if (useLAESA)
                        {
                            laesa = (LAESA<O, T>*) leafnode->index.get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += leafnode->index->getPrunning();
                        }
                        else
                        {
                            for (size_t i = 0; i < leafnode->dataset->getCardinality(); i++)
                            {
                                elementQueue.push(query::ResultEntry<O>(leafnode->dataset->getElement(i).getOID(), this->distanceFunction->operator()(query, leafnode->dataset->getElement(i))));
                            }
                        }

                    }
                    else
                    {
                        innernode = dynamic_cast<btree_type::Inner*>(currentNode);

                        for(size_t i = 0; i < (size_t)(innernode->slotuse + 1); i++)
                        {

                            nodeQueue.push(query::Partition<btree_type::Node*>(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

                        }

                    }

                    if (currentNode->memory_status == MEMORY_STATUS::IN_DISK)
                    {
                        currentNode->clear();
                    }

                }
                else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance())
                {
                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->memory_status == MEMORY_STATUS::IN_DISK)
                    {
                        std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->nodeID);
                        currentNode->deserialize(std::move(nodeData));
                    }

                    if(currentNode->isleafnode()) {

                        leafnode = dynamic_cast<btree_type::Leaf *>(currentNode);
                        this->leafNodeAccess++;

                        if (useLAESA)
                        {
                            laesa = (LAESA<O, T>*) leafnode->index.get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += leafnode->index->getPrunning();
                        }
                        else
                        {
                            for (size_t i = 0; i < leafnode->dataset->getCardinality(); i++)
                            {
                                elementQueue.push(query::ResultEntry<O>(leafnode->dataset->getElement(i).getOID(), this->distanceFunction->operator()(query, leafnode->dataset->getElement(i))));
                            }
                        }

                    }
                    else
                    {
                        innernode = dynamic_cast<btree_type::Inner*>(currentNode);

                        for(size_t i = 0; i < (size_t)(innernode->slotuse + 1); i++)
                        {

                            nodeQueue.push(query::Partition<btree_type::Node*>(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

                        }

                    }

                    if (currentNode->memory_status == MEMORY_STATUS::IN_DISK)
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

        void clear() override
        {
            this->bptree->clear_Tree();
        }

        void print(std::ostream& os) const override
        {

            this->bptree->print(os);

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;
            unsigned long long sz2;

            sz = gervLib::index::Index<O, T>::getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = gervLib::index::Index<O, T>::serialize();
            memcpy(data.get() + offset, indexData.get(), sz);
            offset += sz;
            indexData.reset();

            sz2 = hc->getP();

            memcpy(data.get() + offset, &sz2, sizeof(unsigned long long));
            offset += sizeof(unsigned long long);

            memcpy(data.get() + offset, &numPerLeaf, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numPivots, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numBins, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &storeLeafNode, sizeof(bool));
            offset += sizeof(bool);

            memcpy(data.get() + offset, &storeDirectoryNode, sizeof(bool));
            offset += sizeof(bool);

            memcpy(data.get() + offset, &useLAESA, sizeof(bool));
            offset += sizeof(bool);

            sz = this->keys_serialize->size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, this->keys_serialize->data(), sizeof(unsigned long long) * sz);

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;
            unsigned long long p;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
            memcpy(indexData.get(), _data.get() + offset, sz);
            offset += sz;
            gervLib::index::Index<O, T>::deserialize(std::move(indexData));

            memcpy(&p, _data.get() + offset, sizeof(unsigned long long));
            offset += sizeof(unsigned long long);

            memcpy(&numPerLeaf, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numPivots, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numBins, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&storeLeafNode, _data.get() + offset, sizeof(bool));
            offset += sizeof(bool);

            memcpy(&storeDirectoryNode, _data.get() + offset, sizeof(bool));
            offset += sizeof(bool);

            memcpy(&useLAESA, _data.get() + offset, sizeof(bool));
            offset += sizeof(bool);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            this->keys_serialize = std::make_unique<std::vector<unsigned long long>>(sz);
            memcpy(this->keys_serialize->data(), _data.get() + offset, sizeof(unsigned long long) * sz);
            offset += sizeof(unsigned long long) * sz;

            if (ed != nullptr)
                ed.reset();

            if (hc != nullptr)
                hc.reset();

            if (bptree != nullptr)
            {
                bptree->clear();
                bptree.reset();
            }

            hc = std::make_unique<hilbert::HilbertCurve>(p, numPivots);
            ed = std::make_unique<equidepth::EquiDepth<double>>(numBins, numPivots);
            ed->readFromFile(this->indexFolder);
            bptree = std::make_unique<btree_type>();

            bptree->setHilbertCurve(hc.get());
            bptree->setNumberOfPivots(numPivots);
            bptree->setPageManager(this->pageManager.get());
            bptree->setPivots(this->pivots);
            bptree->setIndexFolder(this->indexFolder);
            bptree->setDistanceFunction(this->distanceFunction.get());

            bptree->buildTreeByLeafKeys(*this->keys_serialize);

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            ans += sizeof(size_t) + gervLib::index::Index<O, T>::getSerializedSize();
            ans += sizeof(unsigned long long);
            ans += sizeof(size_t) * 3;
            ans += sizeof(bool) * 3;

            std::vector<unsigned long long> keys = bptree->getLeafKeys();
            ans += sizeof(size_t);
            ans += sizeof(unsigned long long) * keys.size();

            this->keys_serialize = std::make_unique<std::vector<unsigned long long>>(keys);
            keys.clear();

            return ans;

        }

    private:
        template <typename U>
        bool isInterval(U infBound, U supBound, U test)
        {

            return ((test >= infBound) && (test <= supBound));

        }

        double minDist(std::vector<double>& sq_, unsigned long long** mbr_)
        {

            auto** mbr = new double*[numPivots];

            for(size_t i = 0; i < numPivots; i++)
            {

                mbr[i] = new double[2];
                //mbr[i][0] = numeric_limits<double>::max();
                //mbr[i][1] = numeric_limits<double>::lowest();

                std::pair<double, double> pairMin = ed->getInterval(i, mbr_[i][0]);
                std::pair<double, double> pairMax = ed->getInterval(i, mbr_[i][1]);

                //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
                //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMax.first});
                mbr[i][0] = pairMin.first;

                //mbr[i][1] = std::max({mbr[i][1], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
                //mbr[i][1] = std::max({mbr[i][1], pairMin.second, pairMax.second});
                mbr[i][1] = pairMax.second;


            }

            double limInfCase3 = -1.0;
            double limInfCase2 = -1.0;
            double answer = -1.0;
            bool within = true;

            for(size_t i = 0; i < numPivots; i++)
            {

                if(!isInterval(mbr[i][0], mbr[i][1], sq_[i]))
                {

                    within = false;

                    limInfCase3 = std::max(limInfCase3,
                                           std::min(
                                                   std::abs(sq_[i] - mbr[i][0]),
                                                   std::abs(sq_[i] - mbr[i][1])
                                           )
                    );

                }
                else
                {

                    limInfCase2 = std::min(limInfCase2,
                                           std::min(
                                                   std::abs(sq_[i] - mbr[i][0]),
                                                   std::abs(sq_[i] - mbr[i][1])
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

                if(limInfCase2 != -1)
                {

                    answer = limInfCase2;

                }
                else
                {

                    answer = limInfCase3;

                }

            }

            for(size_t i = 0; i < numPivots; i++)
            {

                delete [] mbr[i];

            }
            delete [] mbr;
//            sq_.clear();

//        cout << "MIN D = " << answer << endl;

            return answer;

        }

        double maxDist(std::vector<double>& sq_, unsigned long long** mbr_)
        {

            auto** mbr = new double*[numPivots];

            for(size_t i = 0; i < numPivots; i++)
            {

                mbr[i] = new double[2];
                //mbr[i][0] = numeric_limits<double>::max();
                //mbr[i][1] = numeric_limits<double>::lowest();

                std::pair<double, double> pairMin = ed->getInterval(i, mbr_[i][0]);
                std::pair<double, double> pairMax = ed->getInterval(i, mbr_[i][1]);

                //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
                //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMax.first});
                mbr[i][0] = pairMin.first;

                //mbr[i][1] = std::max({mbr[i][1], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
                //mbr[i][1] = std::max({mbr[i][1], pairMin.second, pairMax.second});
                mbr[i][1] = pairMax.second;

            }

            double answer = std::numeric_limits<double>::max();

            for(size_t i = 0; i < numPivots; i++)
            {

                if(std::numeric_limits<double>::max() - mbr[i][0] >= sq_[i])
                {

                    answer = std::min(answer, sq_[i] + mbr[i][0]);

                }

                if(std::numeric_limits<double>::max() - mbr[i][1] >= sq_[i])
                {

                    answer = std::min(answer, sq_[i] + mbr[i][1]);

                }


            }

            for(size_t i = 0; i < numPivots; i++)
            {

                delete [] mbr[i];

            }
            delete [] mbr;
//            sq_.clear();

            //cout << "MAX D = " << answer << endl;

            return answer;

        }


    };


}

#endif //GERVLIB_SPBTREE_H
