//
// Created by joaoleite on 9/1/23.
//

#ifndef GERVLIB_SPBTREE_H
#define GERVLIB_SPBTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "btree_multimap.hpp"
#include "EquiDepth.h"
#include "HilbertCurve.h"

namespace gervLib::index::spbtree
{

    template <typename O, typename T, typename H>
    class SPBTree : public Index<O, T>
    {

    public:
        typedef tlx::btree_multimap<H, std::unique_ptr<dataset::BasicArrayObject<O, T>>, O, T, std::less<H>, tlx::btree_default_traits<H, size_t>> btree_type;
        typedef btree_type::Node_t node_type;
        typedef btree_type::LeafNode_t leaf_node_type;
        typedef btree_type::InnerNode_t inner_node_type;

    private:
        std::unique_ptr<equidepth::EquiDepth> equiDepth;
        std::unique_ptr<hilbert::HilbertCurve<H>> hilbertCurve;
        size_t numPerLeaf{}, numPivots{}, num_bins{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::unique_ptr<btree_type> bptree;

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
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->indexFolder = "";
            this->bptree = nullptr;
            this->equiDepth = nullptr;
            this->hilbertCurve = nullptr;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->num_bins = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
        }

        SPBTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
                std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
                std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _num_bins, size_t _numPerLeaf, size_t _pageSize = 0,
                bool _storeDirectoryNode = false, bool _storeLeafNode = false, bool _useLAESA = true, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->numPivots = _numPivots;
            this->numPerLeaf = _numPerLeaf;
            this->pageSize = _pageSize;
            this->num_bins = _num_bins;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->storeLeafNode = _storeLeafNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->bptree = nullptr;

            if (numPivots > 7)
                if constexpr (!configure::is_mpz_class_v<H>)
                    throw std::runtime_error("SPBTree: Number of pivots must be less than 8");

            unsigned long long p = log2(num_bins-1)+1;

            if constexpr (!configure::is_mpz_class_v<H>)
            {
                this->hilbertCurve = std::make_unique<hilbert::HilbertCurve_ull>(p, this->numPivots);
            }
            else
            {
                this->hilbertCurve = std::make_unique<hilbert::HilbertCurve_mpz>(mpz_class(std::to_string(p)), mpz_class(std::to_string(this->numPivots)));
            }

            this->equiDepth = std::make_unique<equidepth::EquiDepth>(this->num_bins, this->numPivots);

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
            this->equiDepth = nullptr;
            this->hilbertCurve = nullptr;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->num_bins = 0;
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
            if (bptree != nullptr) {
                bptree->clear();
                bptree.reset();
            }

            if (equiDepth != nullptr) {
                equiDepth.reset();
            }

            if (hilbertCurve != nullptr) {
                hilbertCurve.reset();
            }

        }

        void buildIndex() override {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);
            this->bptree = std::make_unique<btree_type>();

            std::vector<std::vector<double>> pivotMapping = std::vector<std::vector<double>>(this->dataset->getCardinality(), std::vector<double>(this->numPivots, 0.0));
            std::vector<std::vector<unsigned long long>> discretePivotMapping = std::vector<std::vector<unsigned long long>>(this->dataset->getCardinality(), std::vector<unsigned long long>(this->numPivots, 0));
            std::vector<H> keys(this->dataset->getCardinality());

            std::ofstream file_pivot_mapping(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_pivot_mapping.txt");
            std::ofstream file_disc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_disc.txt");
            std::ofstream file_sfc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_sfc.txt");

            for(size_t i = 0; i < this->dataset->getCardinality(); i++) {

                for (size_t j = 0; j < numPivots; j++) {

                    pivotMapping[i][j] = this->distanceFunction->getDistance(this->dataset->getElement(i), this->pivots->getPivot(j));
                    file_pivot_mapping << pivotMapping[i][j] << " ";

                }

                file_pivot_mapping << "\n";

            }

            this->equiDepth->build(pivotMapping);
            this->equiDepth->saveToFile(this->indexFolder);
            file_pivot_mapping.close();

            for(size_t i = 0; i < this->dataset->getCardinality(); i++)
            {

                if constexpr (configure::is_mpz_class_v<H>)
                {
                    std::vector<mpz_class> tmp = std::vector<mpz_class>(this->numPivots);

                    for(size_t j = 0; j < numPivots; j++)
                    {
                        tmp[j] = mpz_class(this->equiDepth->getBin(j, pivotMapping[i][j]));
                        file_disc << tmp[j] << " ";
                    }

                    file_disc << "\n";

                    keys[i] = this->hilbertCurve->distance_from_point(tmp);
                    file_sfc << keys[i] << "\n";

                    tmp.clear();

                }
                else
                {
                    for(size_t j = 0; j < numPivots; j++)
                    {

                        discretePivotMapping[i][j] = this->equiDepth->getBin(j, pivotMapping[i][j]);
                        file_disc << discretePivotMapping[i][j] << " ";

                    }

                    file_disc << "\n";

                    keys[i] = this->hilbertCurve->distance_from_point(discretePivotMapping[i]);
                    file_sfc << keys[i] << "\n";

                }

            }

            pivotMapping.clear();
            discretePivotMapping.clear();
            this->equiDepth->clear();
            file_disc.close();
            file_sfc.close();

            std::vector<std::pair<H, std::unique_ptr<dataset::BasicArrayObject<O, T>>>> insertValues(keys.size());

            for(size_t i = 0; i < keys.size(); i++)
            {

                std::unique_ptr<dataset::BasicArrayObject<O, T>> obj = std::make_unique<dataset::BasicArrayObject<O, T>>(this->dataset->getElement(i));
                insertValues[i] = std::make_pair(keys[i], std::move(obj));

            }

            keys.clear();

            std::sort(insertValues.begin(), insertValues.end());
            bptree->bulk_load(insertValues.begin(), insertValues.end());

            for (size_t i = 0; i < insertValues.size(); i++) {
                insertValues[i].second->clear();
                insertValues[i].second.reset();
            }

            insertValues.clear();

            bptree->setIndex(this, hilbertCurve.get(), numPivots);
            bptree->build_MBR();

            test();

        }

        void test()
        {
            node_type *root = bptree->getRoot();
        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("Not implemented yet");
        }

        template <typename U>
        bool isInterval(U infBound, U supBound, U test)
        {

            return ((test >= infBound) && (test <= supBound));

        }

        double minDist(dataset::BasicArrayObject<O, T>& auxQuery, std::unique_ptr<std::vector<std::vector<unsigned long long>>>& mbr_)
        {

            std::vector<std::vector<double>> mbr = std::vector<std::vector<double>>(mbr_->size(), std::vector<double>(2, 0.0));

            for (size_t i = 0; i < numPivots; i++)
            {
                std::pair<double, double> pairMin = equiDepth->getInterval(i, (long)mbr_->at(i)[0]);
                std::pair<double, double> pairMax = equiDepth->getInterval(i, (long)mbr_->at(i)[1]);

                mbr[i][0] = pairMin.first;
                mbr[i][1] = pairMax.second;

            }

            double limInfCase3 = -1.0;
            double limInfCase2 = -1.0;
            double answer = -1.0;
            bool within = true;

            for(size_t i = 0; i < numPivots; i++)
            {

                if(!isInterval(mbr[i][0], mbr[i][1], auxQuery[i]))
                {

                    within = false;

                    limInfCase3 = std::max(limInfCase3,
                                           std::min(
                                                   std::abs(auxQuery[i] - mbr[i][0]),
                                                   std::abs(auxQuery[i] - mbr[i][1])
                                           )
                    );

                }
                else
                {

                    limInfCase2 = std::min(limInfCase2,
                                           std::min(
                                                   std::abs(auxQuery[i] - mbr[i][0]),
                                                   std::abs(auxQuery[i] - mbr[i][1])
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

            for (size_t i = 0; i < numPivots; i++)
            {
                mbr[i].clear();
            }

            mbr.clear();

            return answer;

        }

        double maxDist(dataset::BasicArrayObject<O, T>& auxQuery, std::unique_ptr<std::vector<std::vector<unsigned long long>>>& mbr_)
        {

            std::vector<std::vector<double>> mbr = std::vector<std::vector<double>>(mbr_->size(), std::vector<double>(2, 0.0));

            for (size_t i = 0; i < numPivots; i++)
            {
                std::pair<double, double> pairMin = equiDepth->getInterval(i, (long)mbr_->at(i)[0]);
                std::pair<double, double> pairMax = equiDepth->getInterval(i, (long)mbr_->at(i)[1]);

                mbr[i][0] = pairMin.first;
                mbr[i][1] = pairMax.second;

            }

            double answer = std::numeric_limits<double>::max();

            for(size_t i = 0; i < numPivots; i++)
            {

                if(std::numeric_limits<double>::max() - mbr[i][0] >= auxQuery[i])
                {

                    answer = std::min(answer, auxQuery[i] + mbr[i][0]);

                }

                if(std::numeric_limits<double>::max() - mbr[i][1] >= auxQuery[i])
                {

                    answer = std::min(answer, auxQuery[i] + mbr[i][1]);

                }


            }

            for (size_t i = 0; i < numPivots; i++)
            {
                mbr[i].clear();
            }

            mbr.clear();

            return answer;

        }


        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            this->leafNodeAccess = 0;
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            std::priority_queue<query::Partition<node_type*>, std::vector<query::Partition<node_type*>>, std::greater<query::Partition<node_type*>>> nodeQueue;
            std::priority_queue<query::ResultEntry<O>, std::vector<query::ResultEntry<O>>, std::greater<query::ResultEntry<O>>> elementQueue;
            query::Result<O> result;
            result.setMaxSize(k);
            query::Partition<node_type*> currentPartition;
            node_type* currentNode;
            leaf_node_type* currentLeafNode;
            inner_node_type* currentInnerNode;

            nodeQueue.push(query::Partition<node_type*>(bptree->getRoot(), 0.0, std::numeric_limits<double>::max()));

            dataset::BasicArrayObject<O, T> auxQuery = query;

            for (size_t i = 0; i < numPivots; i++)
            {
                auxQuery.operator[](i) = this->distanceFunction->getDistance(query, this->pivots->getPivot(i));
            }

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty()))
            {

                if (elementQueue.empty()) {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->is_leafnode())
                    {

                        currentLeafNode = static_cast<leaf_node_type*>(currentNode);

                        for (size_t i = 0; i < currentNode->slotuse; i++)
                        {
                            elementQueue.push(query::ResultEntry<O>(currentLeafNode->slotdata[i].second->getOID(), this->distanceFunction->getDistance(query, *currentLeafNode->slotdata[i].second)));
                        }

                    }
                    else
                    {

                        currentInnerNode = static_cast<inner_node_type*>(currentNode);

                        for (size_t i = 0; i < currentNode->slotuse + 1; i++)
                        {

                            nodeQueue.push(query::Partition<node_type*>(currentInnerNode->childid[i], minDist(auxQuery, currentInnerNode->mbr), maxDist(auxQuery, currentInnerNode->mbr)));

                        }

                    }

                }
                else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance())
                {
                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->is_leafnode())
                    {

                        currentLeafNode = static_cast<leaf_node_type*>(currentNode);

                        for (size_t i = 0; i < currentNode->slotuse; i++)
                        {
                            elementQueue.push(query::ResultEntry<O>(currentLeafNode->slotdata[i].second->getOID(), this->distanceFunction->getDistance(query, *currentLeafNode->slotdata[i].second)));
                        }

                    }
                    else
                    {

                        currentInnerNode = static_cast<inner_node_type*>(currentNode);

                        for (size_t i = 0; i < currentNode->slotuse + 1; i++)
                        {

                            nodeQueue.push(query::Partition<node_type*>(currentInnerNode->childid[i], minDist(auxQuery, currentInnerNode->mbr), maxDist(auxQuery, currentInnerNode->mbr)));

                        }

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


    };

}

#endif //GERVLIB_SPBTREE_H
