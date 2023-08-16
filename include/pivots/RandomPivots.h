//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_RANDOMPIVOTS_H
#define GERVLIB_RANDOMPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{


    template <typename O, typename T>
    class RandomPivots : public Pivot<O, T>
    {

    public:
        RandomPivots() : Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::RANDOM); }
        RandomPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset, size_t nPivots, size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::RANDOM);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, nullptr, nPivots);
        }

        ~RandomPivots() override = default;

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
                    throw std::runtime_error("RandomPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
            }
            else
            {

                utils::Timer timer;
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

                std::vector<size_t> indexes = utils::generateRandomNumbers(0, sample->getCardinality()-1, nPivots, false, this->getSeed());

                for (size_t i = 0; i < nPivots; i++)
                {

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << i << ": " << sample->operator[](indexes[i]) << std::endl;
#endif

                    this->setPivot(i, sample->getElement(indexes[i]));

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                indexes.clear();

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

#endif //GERVLIB_RANDOMPIVOTS_H
