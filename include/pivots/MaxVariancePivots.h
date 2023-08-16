//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_MAXVARIANCEPIVOTS_H
#define GERVLIB_MAXVARIANCEPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class MaxVariancePivots: public Pivot<O, T>
    {

    public:

        MaxVariancePivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::MAXVARIANCE); }

        MaxVariancePivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                          std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                          size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::MAXVARIANCE);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~MaxVariancePivots() override = default;

        void operator()(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots) override
        {

            if (nPivots >= dataset->getCardinality())
            {
                if (nPivots == dataset->getCardinality())
                {
                    this->pivots->setData(dataset->getData());
                    return;
                }
                else
                    throw std::runtime_error("MaxVariancePivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
            }
            else
            {

                utils::Timer timer{};
                timer.start();

                nPivots += this->getNumberOfDropPivots();
                this->setNumberOfPivots(nPivots);
                this->pivots->setPath(dataset->getPath());

                bool smpl = this->needSample(dataset);
                std::unique_ptr<dataset::Dataset<O, T>> sample;

                if (smpl)
                    sample = dataset->sample(this->getSampleSize(), false, this->getSeed());
                else
                    sample = std::move(dataset);

                std::vector<size_t> aux(sample->getCardinality());
                utils::Random<size_t> rGen(this->getSeed());
                std::iota(aux.begin(), aux.end(), 0);
                std::shuffle(aux.begin(), aux.end(), rGen.generator);
                std::vector<size_t> cand(aux.begin(), aux.begin() + aux.size()/2),
                        elements(aux.begin() + aux.size()/2, aux.end());

                double mean, sum, dist;
                std::vector<std::pair<double, dataset::BasicArrayObject<O, T>>> maxvar;

#ifdef ENABLE_DEBUG

                std::cout << "====================================================================================\n";

                for (size_t x = 0; x < cand.size(); x++)
                    std::cout << "Candidate " << x << ": " << sample->operator[](cand[x]) << std::endl;

                std::cout << "\n====================================================================================\n\n";

                for (size_t x = 0; x < elements.size(); x++)
                    std::cout << "Element " << x << ": " << sample->operator[](elements[x]) << std::endl;

                std::cout << "====================================================================================\n\n";

#endif

                for(size_t x = 0; x < cand.size(); x++)
                {

                    mean = 0.0, sum = 0.0;

                    for(size_t y = 0; y < elements.size(); y++)
                        mean += df->operator()(sample->getElement(cand[x]), sample->getElement(elements[y]));

                    mean /= static_cast<double>(elements.size());

                    for(size_t y = 0; y < elements.size(); y++)
                    {

                        dist = df->operator()(sample->getElement(cand[x]), sample->getElement(elements[y]));
                        sum += (dist - mean) * (dist - mean);

                    }

                    sum /= static_cast<double>(elements.size());
                    maxvar.push_back(std::make_pair(sum, sample->getElement(cand[x])));

                }

                std::sort(maxvar.begin(), maxvar.end(), std::greater());

                for(size_t x = 0; x < nPivots; x++)
                {

                    this->setPivot(x, maxvar[x].second);

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << x << ": " << maxvar[x].second << std::endl;
#endif

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                aux.clear();
                cand.clear();
                elements.clear();
                maxvar.clear();

                if (this->drop != 0)
                {
                    for (size_t i = 0; i < this->drop; i++)
                        this->pivots->erase(0);
                }

                timer.stop();
                this->incrementElapsedTime(timer.getElapsedTime());

            }

        }

    };

}

#endif //GERVLIB_MAXVARIANCEPIVOTS_H
