//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_SSSPIVOTS_H
#define GERVLIB_SSSPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class SSSPivots: public Pivot<O, T>
    {

    private:
        inline static double alpha = 0.0015;

    public:
        SSSPivots() : Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::SSS); }

        SSSPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                  std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                  size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::SSS);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~SSSPivots() override = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            memcpy(data.get() + offset, &this->alpha, sizeof(double));
            offset += sizeof(double);

            memcpy(data.get() + offset, &this->seed, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &this->drop, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &this->elapsedTime, sizeof(long long));
            offset += sizeof(long long);

            size = this->sampleSize.index();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (this->sampleSize.index() == 0)
            {

                memcpy(data.get() + offset, &std::get<size_t>(this->sampleSize), sizeof(size_t));
                offset += sizeof(size_t);

            }
            else
            {

                memcpy(data.get() + offset, &std::get<double>(this->sampleSize), sizeof(double));
                offset += sizeof(double);

            }

            size = this->pivots == nullptr ? 0 : this->pivots->getSerializedSize();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (this->pivots != nullptr)
            {

                std::unique_ptr<u_char[]> aux = this->pivots->serialize();
                memcpy(data.get() + offset, aux.get(), size);
                offset += size;
                aux.reset();

            }

            size = PIVOT_TYPE2STR[this->type].size();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, PIVOT_TYPE2STR[this->type].c_str(), size);
            offset += size;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, size;

            memcpy(&this->alpha, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            memcpy(&this->seed, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&this->drop, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&this->elapsedTime, _data.get() + offset, sizeof(long long));
            offset += sizeof(long long);

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                size_t sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);
                this->sampleSize = sz2;

            }
            else
            {

                double sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                this->sampleSize = sz2;

            }

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                this->clear();

            }
            else
            {

                this->pivots->clear();
                this->pivots = std::make_unique<dataset::Dataset<O, T>>();
                std::unique_ptr<u_char[]> aux(new u_char[size]);
                memcpy(aux.get(), _data.get() + offset, size);
                offset += size;
                this->pivots->deserialize(std::move(aux));

            }

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string str;
            str.resize(size);
            memcpy(&str[0], _data.get() + offset, size);
            offset += size;

            _data.reset();

        }

        size_t getSerializedSize() override
        {

            size_t ans = sizeof(double) + //alpha
                         sizeof(size_t) + //seed
                         sizeof(size_t) + //drop
                         sizeof(long long) + //time
                         sizeof(size_t) + (this->sampleSize.index() == 0 ? sizeof(size_t) : sizeof(double)); //number to indicate variant possibility and variable size

            ans += sizeof(size_t) + (this->pivots == nullptr ? 0 : this->pivots->getSerializedSize());
            std::string aux = PIVOT_TYPE2STR[this->getPivotType()];
            ans += sizeof(size_t) + aux.size();

            return ans;

        }

        void print(std::ostream& os) const override
        {

            os << "Pivot Type: " << this->type << std::endl;
            os << "Pivot Seed: " << this->seed << std::endl;
            os << "Pivot Elapsed Time: " << this->elapsedTime << std::endl;
            os << "Alpha: " << this->alpha << std::endl;
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Number of drop pivots: " << this->drop << std::endl;
            os << "Pivot Number of Pivots: " << this->pivots->getCardinality() << std::endl;
            os << "Pivot Pivots: " << std::endl;
            os << *this->pivots << std::endl;

        }

        double getAlpha(){ return alpha; }

        void setAlpha(double _alpha){ alpha = _alpha; }

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
                    throw std::runtime_error("SSSPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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

                size_t currentPivot = 0, count = 0, p1;
                std::vector<bool> bitmap(sample->getCardinality(), false);
                std::vector<size_t> pvtIndex(nPivots, 0);
                std::vector<size_t> aux = utils::generateRandomNumbers(0, sample->getCardinality()-1, 1, false, this->getSeed());
                double max = std::numeric_limits<double>::min(), dist, threshold;

                p1 = aux[0];
                bitmap[p1] = true;
                this->setPivot(currentPivot, sample->getElement(p1));

#ifdef ENABLE_DEBUG
                std::cout << "Pivot " << currentPivot << ": " << sample->operator[](p1) << std::endl;
#endif

                pvtIndex[currentPivot++] = p1;

                for(size_t x = 0; x < sample->getCardinality(); x++)
                {

                    for(size_t y = x + 1; y < sample->getCardinality(); y++)
                    {

                        dist = df->operator()(sample->getElement(x), sample->getElement(y));
                        max = std::max(max, dist);

                    }

                }

                threshold = alpha * max;

                for(size_t x = 0; x < sample->getCardinality(); x++)
                {

                    if (!bitmap[x])
                    {

                        count = 0;

                        for(size_t y = 0; y < currentPivot; y++)
                        {

                            dist = df->operator()(sample->getElement(x), sample->getElement(pvtIndex[y]));

                            if(dist >= threshold)
                                count++;

                        }

                        if (count == currentPivot)
                        {

                            bitmap[x] = true;
                            this->setPivot(currentPivot, sample->getElement(x));

#ifdef ENABLE_DEBUG
                            std::cout << "Pivot " << currentPivot << ": " << sample->operator[](x) << std::endl;
#endif

                            pvtIndex[currentPivot++] = x;

                        }

                        if(currentPivot == nPivots)
                            break;

                    }

                }

                //Fill with random pivots
                if (currentPivot < nPivots)
                {

                    utils::Random<size_t> rGen(this->getSeed());
                    size_t pos;

                    while (currentPivot < nPivots)
                    {

                        pos = rGen(0, sample->getCardinality());

                        if(!bitmap[pos]) {

                            if(pvtIndex.end() == std::find(pvtIndex.begin(), pvtIndex.end(), pos))
                            {

                                bitmap[pos] = true;
                                this->setPivot(currentPivot, sample->getElement(pos));

#ifdef ENABLE_DEBUG
                                std::cout << "Random Pivot " << currentPivot << ": " << sample->operator[](pos) << std::endl;
#endif

                                pvtIndex[currentPivot++] = pos;

                            }

                        }

                    }

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

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


#endif //GERVLIB_SSSPIVOTS_H
