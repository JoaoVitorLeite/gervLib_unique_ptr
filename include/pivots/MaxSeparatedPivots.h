//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_MAXSEPARATEDPIVOTS_H
#define GERVLIB_MAXSEPARATEDPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class MaxSeparatedPivots : public Pivot<O, T>
    {

    public:
        MaxSeparatedPivots() : Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::MAXSEPARATED); }
        MaxSeparatedPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                           std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                           size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::MAXSEPARATED);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~MaxSeparatedPivots() override = default;

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
                    throw std::runtime_error("MaxSeparatedPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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

                size_t currentPivot = 0, pos = 0, p1;
                std::vector<bool> bitmap(sample->getCardinality(), false);
                std::vector<size_t> pvtIndex(nPivots, 0);
                std::vector<size_t> aux = utils::generateRandomNumbers(0, sample->getCardinality()-1, 1, false, this->getSeed());
                double max = std::numeric_limits<double>::min(), sum = 0.0, dist;

                p1 = aux[0];
                bitmap[p1] = true;
                this->setPivot(currentPivot, sample->getElement(p1));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << currentPivot << ": " << sample->operator[](p1) << std::endl;
#endif

                pvtIndex[currentPivot++] = p1;


                if (nPivots == 1)
                {
                    this->pivots->setDimensionality(sample->getDimensionality());

                    if (smpl)
                        sample.reset();
                    else
                        dataset = std::move(sample);

                    bitmap.clear();
                    pvtIndex.clear();
                    aux.clear();

                    if (this->drop != 0)
                    {
                        for (size_t i = 0; i < this->drop; i++)
                            this->pivots->erase(0);
                    }

                    timer.stop();
                    this->incrementElapsedTime(timer.getElapsedTime());

                    return;
                }

                for (size_t x = 0; x < sample->getCardinality(); x++)
                {

                    dist = df->operator()(sample->getElement(x), sample->getElement(p1));

                    if (dist > max)
                    {

                        max = dist;
                        pos = x;

                    }

                }

                bitmap[pos] = true;
                this->setPivot(currentPivot, sample->getElement(pos));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << currentPivot << ": " << sample->operator[](pos) << std::endl;
#endif

                pvtIndex[currentPivot++] = pos;

                if (nPivots == 2)
                {
                    this->pivots->setDimensionality(sample->getDimensionality());

                    if (smpl)
                        sample.reset();
                    else
                        dataset = std::move(sample);

                    bitmap.clear();
                    pvtIndex.clear();
                    aux.clear();

                    if (this->drop != 0)
                    {
                        for (size_t i = 0; i < this->drop; i++)
                            this->pivots->erase(0);
                    }

                    timer.stop();
                    this->incrementElapsedTime(timer.getElapsedTime());

                    return;
                }

                while (currentPivot < nPivots)
                {

                    max = std::numeric_limits<double>::min();

                    for (size_t x = 0; x < sample->getCardinality(); x++)
                    {

                        if (!bitmap[x])
                        {

                            sum = 0.0;

                            for(size_t y = 0; y < currentPivot; y++)
                                sum += df->operator()(sample->getElement(x), sample->getElement(pvtIndex[y]));

                            if (sum > max)
                            {

                                max = sum;
                                pos = x;

                            }

                        }

                    }

                    bitmap[pos] = true;
                    this->setPivot(currentPivot, sample->getElement(pos));

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << currentPivot << ": " << sample->operator[](pos) << std::endl;
#endif

                    pvtIndex[currentPivot++] = pos;

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                bitmap.clear();
                pvtIndex.clear();
                aux.clear();

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

#endif //GERVLIB_MAXSEPARATEDPIVOTS_H
