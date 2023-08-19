//
// Created by joaovictor on 19/08/23.
//

#ifndef GERVLIB_SEQUENTIALSCAN_H
#define GERVLIB_SEQUENTIALSCAN_H

#include "Index.h"

namespace gervLib::index
{

    template <typename O, typename T>
    class SequentialScan : public Index<O, T>
    {

    protected:
        std::string headerBuildFile() override
        {
            return "";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount";
        }

    public:
        SequentialScan()
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::SEQUENTIAL_SCAN;
            this->indexName = "SEQUENTIAL_SCAN";
            this->indexFolder = "";
        }

        SequentialScan(std::unique_ptr<dataset::Dataset<O, T>> _dataset, std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df, std::string folder="")
        {
            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::SEQUENTIAL_SCAN;
            this->indexName = "SEQUENTIAL_SCAN";

            if (!folder.empty())
                this->indexFolder = folder;

            this->buildIndex();
            this->generateIndexFiles(true, true);

        }

        explicit SequentialScan(std::string _folder, std::string serializeFile = "")
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::SEQUENTIAL_SCAN;
            this->indexName = "SEQUENTIAL_SCAN";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (!serializeFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializeFile);
        }

        ~SequentialScan() override = default;

        void buildIndex() override { }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            return this->indexType == other->getIndexType() &&
                   this->dataset->isEqual(*other->getDataset()) &&
                   this->distanceFunction->isEqual(other->getDistanceFunction());
        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();

            std::vector<gervLib::query::ResultEntry<O>> results;

            for (size_t i = 0; i < this->dataset->getCardinality(); i++)
            {
                double dist = this->distanceFunction->getDistance(query, this->dataset->getElement(i));
                results.emplace_back(i, dist);
            }

            std::sort(results.begin(), results.end(),
                      [](const gervLib::query::ResultEntry<O> &a, const gervLib::query::ResultEntry<O> &b) {
                          return a.getDistance() < b.getDistance();
                      });
            results.resize(k);

            std::string expt_id = utils::generateExperimentID();
            timer.stop();

            if (saveResults)
            {

                this->saveResultToFile(results, query, "kNN", expt_id, {expt_id, std::to_string(k), "-1",
                                                                        std::to_string(timer.getElapsedTime()),
                                                                        std::to_string(timer.getElapsedTimeSystem()),
                                                                        std::to_string(timer.getElapsedTimeUser()),
                                                                        std::to_string(this->distanceFunction->getDistanceCount())});


            }

            return results;

        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("SequentialScan::kNNIncremental not implemented yet");
        }


    };

}

#endif //GERVLIB_SEQUENTIALSCAN_H
