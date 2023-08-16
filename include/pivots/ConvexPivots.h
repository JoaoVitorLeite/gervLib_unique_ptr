//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_CONVEXPIVOTS_H
#define GERVLIB_CONVEXPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class ConvexPivots: public Pivot<O, T>
    {

    public:
        ConvexPivots() : Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::CONVEX); }
        ConvexPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                     std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                     size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::CONVEX);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~ConvexPivots() override = default;

        virtual void operator()(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
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
                double max = std::numeric_limits<double>::min(), dist;
                std::vector<bool> bitmap(sample->getCardinality(), false);
                std::vector<size_t> pvtIndex(nPivots, 0);
                std::vector<size_t> aux = utils::generateRandomNumbers(0, sample->getCardinality()-1, nPivots, false,
                                                                       this->getSeed());
                p1 = aux[0];
                bitmap[p1] = true;
                pvtIndex[currentPivot++] = p1;
                this->setPivot(0, sample->getElement(p1));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << 0 << ": " << sample->operator[](p1) << std::endl;
#endif

                if(nPivots == 1)
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

                for(size_t x = 0; x < sample->getCardinality(); x++)
                {

                    dist = df->operator()(sample->getElement(p1), sample->getElement(x));

                    if (dist > max)
                    {

                        max = dist;
                        pos = x;

                    }

                }

                bitmap[pos] = true;
                pvtIndex[currentPivot++] = pos;
                this->setPivot(1, sample->getElement(pos));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << 1 << ": " << sample->operator[](pos) << std::endl;
#endif

                if(nPivots == 2)
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

                max = std::numeric_limits<double>::min();
                pos = 0;

                for(size_t x = 0; x < sample->getCardinality(); x++)
                {

                    dist = df->operator()(sample->getElement(pvtIndex[1]), sample->getElement(x));

                    if ((dist > max) && (!bitmap[x]))
                    {

                        max = dist;
                        pos = x;

                    }

                }

                bitmap[pos] = true;
                pvtIndex[currentPivot++] = pos;
                this->setPivot(2, sample->getElement(pos));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << 2 << ": " << sample->operator[](pos) << std::endl;
#endif

                if(nPivots == 3)
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

                double edge = df->operator()(sample->getElement(pvtIndex[1]), sample->getElement(pvtIndex[2]));

                while (currentPivot < nPivots)
                {

                    double error = 0;
                    double min = std::numeric_limits<double>::max();
                    pos = 0;

                    for(size_t x = 0; x < sample->getCardinality(); x++)
                    {

                        if(!bitmap[x])
                        {

                            error = 0.0;

                            for(size_t y = 0; y < currentPivot; y++)
                            {

                                error += std::abs(edge - df->operator()(sample->getElement(pvtIndex[y]), sample->getElement(x)));

                            }

                            if(error < min)
                            {

                                min = error;
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

#endif //GERVLIB_CONVEXPIVOTS_H
