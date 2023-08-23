//
// Created by joaovictor on 20/08/23.
//

#ifndef GERVLIB_LAESA_H
#define GERVLIB_LAESA_H

#include "Index.h"
#include "Matrix.h"

namespace gervLib::index
{

    template <typename O, typename T>
    class LAESA : public Index<O, T>
    {

    private:
        std::unique_ptr<matrix::Matrix<O, T>> matrix;

    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning";
        }

    public:
        LAESA()
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->matrix = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::LAESA_t;
            this->indexName = "LAESA";
            this->indexFolder = "";
        }

        LAESA(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
              std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
              std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t numPivots, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->pivots->operator()(this->dataset, this->distanceFunction, numPivots);
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::LAESA_t;
            this->indexName = "LAESA";
            this->matrix = std::make_unique<matrix::Matrix<O, T>>(numPivots, this->dataset->getCardinality());

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);
            this->buildIndex();

        }

        explicit LAESA(std::string _folder, std::string serializedFile = "")
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->matrix = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::LAESA_t;
            this->indexName = "LAESA";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~LAESA() override
        {
            if (matrix != nullptr)
                matrix.reset();
        }

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();

            for(size_t j = 0; j < this->dataset->getCardinality(); j++)
            {

                for(size_t i = 0; i < this->pivots->getNumberOfPivots(); i++)
                {

                    matrix->setValue(i, j, this->distanceFunction->operator()(this->pivots->getPivot(i), this->dataset->getElement(j)));

                }

            }

            timer.stop();

            std::ofstream buildFile(this->buildFile, std::ios::app);

            if (buildFile.is_open())
            {
                buildFile << timer.getElapsedTimeSystem()
                          << ","
                          << timer.getElapsedTimeSystem()
                          << ","
                          << timer.getElapsedTimeUser()
                          << ","
                          << this->distanceFunction->getDistanceCount()
                          << std::endl;
            }
            else
                throw std::runtime_error("LAESA::buildIndex: Error opening file " + this->buildFile);

            buildFile.close();

        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {

            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<LAESA<O, T>*>(other.get());
            std::unique_ptr<LAESA<O, T>> _derived;

            if (_other == nullptr)
                return false;
            else
            {
                other.release();
                _derived.reset(_other);
                if (!matrix->isEqual(_derived->matrix))
                    return false;
            }

            other = std::move(_derived);

            return true;
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("LAESA::kNNIncremental not implemented yet");
        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            std::vector<query::ResultEntry<O>> results;
            std::vector<std::pair<double, O>> shapiro;
            std::vector<bool> bitmap(this->dataset->getCardinality(), false);
            gervLib::query::Result<O> resultPQ(false, k);

            for (size_t i = 0; i < this->pivots->getNumberOfPivots(); i++)
                shapiro.push_back(std::make_pair(this->distanceFunction->operator()(query, this->pivots->getPivot(i)), i));

            std::sort(shapiro.begin(), shapiro.end());

            double dist, distToPivot, lowerBound, upperBound, radius = std::numeric_limits<double>::max();
            size_t p;

            for (size_t x = 0; x < shapiro.size(); x++)
            {

                p = shapiro[x].second;

                if(!bitmap[this->pivots->getPivot(p).getOID()])
                {

                    distToPivot = shapiro[x].first;

                    for (size_t i = 0; i < this->dataset->getCardinality(); i++)
                    {

                        if(!bitmap[i])
                        {

                            lowerBound = std::abs(distToPivot - matrix->getValue(p, i));
                            upperBound = distToPivot + matrix->getValue(p, i);

                            if (upperBound <= radius)
                            {

                                resultPQ.push(gervLib::query::ResultEntry<O>(this->dataset->getElement(i).getOID(), this->distanceFunction->operator()(query, this->dataset->getElement(i))));
                                bitmap[i] = true;

                                if (resultPQ.size() >= k)
                                    radius = resultPQ.top().getDistance();

                            }

                            if (lowerBound > radius)
                            {

                                bitmap[i] = true;
                                this->prunning++;

                            }

                        }

                    }

                }

            }

            for (size_t i = 0; i < this->dataset->getCardinality(); i++)
            {

                if (!bitmap[i])
                {

                    dist = this->distanceFunction->operator()(query, this->dataset->getElement(i));

                    if (dist < radius)
                    {

                        resultPQ.push(gervLib::query::ResultEntry<O>(this->dataset->getElement(i).getOID(), dist));
                        radius = resultPQ.top().getDistance();

                    }

                }

            }

            results = resultPQ.getResults();
            std::reverse(results.begin(), results.end());
            resultPQ.clear();

            std::string expt_id = gervLib::utils::generateExperimentID();

            if (saveResults)
            {

                this->saveResultToFile(results, query, "kNN", expt_id, {expt_id, std::to_string(k), "-1",
                                                                        std::to_string(timer.getElapsedTime()),
                                                                        std::to_string(timer.getElapsedTimeSystem()),
                                                                        std::to_string(timer.getElapsedTimeUser()),
                                                                        std::to_string(this->distanceFunction->getDistanceCount()),
                                                                        std::to_string(this->prunning)});

            }

            return results;

        }

        void print(std::ostream& os) const override
        {
            os << "Index Name: " << this->indexName << std::endl;
            os << "Index Type: " << indexTypeMap.at(this->indexType) << std::endl;
            os << "Index Folder: " << this->indexFolder << std::endl;
            os << "Build File: " << this->buildFile << std::endl;
            os << "Experiment File: " << this->experimentFile << std::endl;
            os << "Page Size: " << this->pageSize << std::endl;

            if (this->dataset != nullptr)
                os << "Dataset: " << std::endl << *this->dataset << std::endl;
            else
                os << "Dataset: " << std::endl << "NULL" << std::endl;

            if (this->distanceFunction != nullptr)
                os << "Distance Function: " << std::endl << *this->distanceFunction << std::endl;
            else
                os << "Distance Function: " << std::endl << "NULL" << std::endl;

            if (this->pivots != nullptr)
                os << "Pivots: " << std::endl << *this->pivots << std::endl;
            else
                os << "Pivots: " << std::endl << "NULL" << std::endl;

            if (this->pageManager != nullptr)
                os << "Page Manager: " << std::endl << *this->pageManager << std::endl;
            else
                os << "Page Manager: " << std::endl << "NULL" << std::endl;

            if (matrix != nullptr)
                os << "Matrix: " << std::endl << *matrix << std::endl;
            else
                os << "Matrix: " << std::endl << "NULL" << std::endl;
        }

        void clear() override
        {
            index::Index<O, T>::clear();

            if (matrix != nullptr)
                matrix.reset();

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, sz;

            if (matrix != nullptr)
            {
                sz = matrix->getSerializedSize();
                std::memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> _data = matrix->serialize();
                std::memcpy(data.get() + offset, _data.get(), sz);
                offset += sz;
                _data.reset();

            }
            else
            {
                sz = 0;
                std::memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            sz = index::Index<O, T>::getSerializedSize();
            std::memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> _data = index::Index<O, T>::serialize();
            std::memcpy(data.get() + offset, _data.get(), sz);
            offset += sz;
            _data.reset();

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            std::memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0)
            {
                matrix = std::make_unique<matrix::Matrix<O, T>>();
                std::unique_ptr<u_char[]> _matrix(new u_char[sz]);
                std::memcpy(_matrix.get(), _data.get() + offset, sz);
                offset += sz;
                matrix->deserialize(std::move(_matrix));
            }
            else
                matrix = nullptr;

            std::memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> _index(new u_char[sz]);
            std::memcpy(_index.get(), _data.get() + offset, sz);
            offset += sz;
            index::Index<O, T>::deserialize(std::move(_index));

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            return sizeof(size_t)*2 +
                   gervLib::index::Index<O, T>::getSerializedSize() +
                   ((matrix != nullptr) ? matrix->getSerializedSize() : 0);
        }

    };


}


#endif //GERVLIB_LAESA_H
