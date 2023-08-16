//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_HFIPIVOTS_H
#define GERVLIB_HFIPIVOTS_H

#include "DatasetWrapper.h"
#include "ConvexPivots.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class HFIPivots : public Pivot<O, T>
    {

    private:
        inline static std::variant<size_t, double> pivotSize = -1.0, limInfSampleSize = (size_t)300;
        inline static size_t convexDropPivots = 0;

    public:
        HFIPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::HFI); }

        HFIPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                  std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                  size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::HFI);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~HFIPivots() override = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            memcpy(data.get() + offset, &this->convexDropPivots, sizeof(size_t));
            offset += sizeof(size_t);

            size = this->pivotSize.index();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (this->pivotSize.index() == 0)
            {

                memcpy(data.get() + offset, &std::get<size_t>(this->pivotSize), sizeof(size_t));
                offset += sizeof(size_t);

            }
            else
            {

                memcpy(data.get() + offset, &std::get<double>(this->pivotSize), sizeof(double));
                offset += sizeof(double);

            }

            size = this->limInfSampleSize.index();

            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (this->limInfSampleSize.index() == 0)
            {

                memcpy(data.get() + offset, &std::get<size_t>(this->limInfSampleSize), sizeof(size_t));
                offset += sizeof(size_t);

            }
            else
            {

                memcpy(data.get() + offset, &std::get<double>(this->limInfSampleSize), sizeof(double));
                offset += sizeof(double);

            }

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

            memcpy(&this->convexDropPivots, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                size_t sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);
                this->pivotSize = sz2;

            }
            else
            {

                double sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                this->pivotSize = sz2;

            }

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                size_t sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);
                this->limInfSampleSize = sz2;

            }
            else
            {

                double sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                this->limInfSampleSize = sz2;

            }

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

            size_t ans = sizeof(size_t) + //drop convex pivots
                         sizeof(size_t) + (this->pivotSize.index() == 0 ? sizeof(size_t) : sizeof(double)) + //pivot size
                         sizeof(size_t) + (this->limInfSampleSize.index() == 0 ? sizeof(size_t) : sizeof(double)) + //lim inf sample size
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
            os << "Number of pairs: " << configure::variant2string(this->limInfSampleSize) << std::endl;
            os << "Pivot candidates size: " << configure::variant2string(this->pivotSize) << std::endl;
            os << "Number of drop convex pivots: " << this->convexDropPivots << std::endl;
            os << "Number of drop pivots: " << this->drop << std::endl;
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Pivot Number of Pivots: " << this->pivots->getCardinality() << std::endl;
            os << "Pivot Pivots: " << std::endl;
            os << *this->pivots << std::endl;

        }

        void setPivotSize(size_t sz){ pivotSize = sz; }

        void setPivotSize(double sz){ pivotSize = sz; }

        std::variant<size_t, double> getPivotSize(){ return pivotSize; }

        size_t getConvexDropPivots(){ return convexDropPivots; }

        void setConvexDropPivots(size_t _drop){ convexDropPivots = _drop; }

        void setLimInfSampleSize(size_t sz){ limInfSampleSize = sz; }

        void setLimInfSampleSize(double sz){ limInfSampleSize = sz; }

        std::variant<size_t, double> getLimInfSampleSize(){ return limInfSampleSize; }

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
                    throw std::runtime_error("HFIPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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

                size_t pvt_sz = 0, limInf = 0;

                if(pivotSize.index() == 0) {
                    pvt_sz = std::get<size_t>(pivotSize);
                }
                else {
                    if (std::get<double>(pivotSize) == -1.0)
                        pvt_sz = nPivots * 2;
                    else
                        pvt_sz = std::floor(std::get<double>(pivotSize) * sample->getCardinality());
                }

                if (pvt_sz <= nPivots)
                    throw std::invalid_argument("HFIPivots::operator(): Pivot size must be greater than the number of pivots.");

                if(limInfSampleSize.index() == 0)
                    limInf = std::get<size_t>(limInfSampleSize);
                else
                    limInf = std::floor(std::get<double>(limInfSampleSize) * sample->getCardinality());

                if (pvt_sz > sample->getCardinality())
                    pvt_sz = sample->getCardinality();

                std::vector<bool> bitmap(pvt_sz, false);
                std::vector<size_t> pivot_index(pvt_sz);
                std::iota(pivot_index.begin(), pivot_index.end(), 0);
                std::vector<std::vector<size_t>> pairs_index(limInf, std::vector<size_t>(2));
                utils::Random<size_t> rGen(this->getSeed());
                std::vector<double> dist_pairs(limInf);
                std::vector<std::vector<double>> dist_pivots_pairs(pvt_sz, std::vector<double>(limInf, 0.0));

                std::unique_ptr<dataset::DatasetWrapper<O, T>> wrapper = std::make_unique<dataset::DatasetWrapper<O, T>>(sample);
                ConvexPivots<O, T> convex = ConvexPivots<O, T>();
                convex.setSeed(this->getSeed());

#ifdef ENABLE_DEBUG
                std::cout << "============= Convex candidates (reset index): " << pvt_sz << " =============================" << std::endl;
#endif

                size_t old_drop = this->getNumberOfDropPivots();
                convex.setNumberOfDropPivots(convexDropPivots);
                convex(sample, df, pvt_sz);
                this->setNumberOfDropPivots(old_drop);

                for(size_t x = 0; x < pvt_sz; x++)
                    pivot_index[x] = convex.getPivot(x).getOID();

                convex.setNumberOfPivots(0);

#ifdef ENABLE_DEBUG
                std::cout << "==============================================================================\n" << std::endl;
                std::cout << "============= Pairs candidates: " << limInf << " ============================================\n" << std::endl;
#endif

                for(size_t x = 0; x < limInf; x++)
                {

                    while (pairs_index[x][0] == pairs_index[x][1])
                    {

                        pairs_index[x][0] = rGen(0, sample->getCardinality()-1);
                        pairs_index[x][1] = rGen(0, sample->getCardinality()-1);

                    }

#ifdef ENABLE_DEBUG
                    std::cout << "Pair " << x << ": " << sample->operator[](pairs_index[x][0]) << " <->  " << sample->operator[](pairs_index[x][1]) << std::endl;
#endif

                    dist_pairs[x] = df->operator()(sample->getElement(pairs_index[x][0]), sample->getElement(pairs_index[x][1]));

                    for(size_t y = 0; y < pvt_sz; y++)
                    {

                        dist_pivots_pairs[y][x] = fabs(df->operator()(sample->getElement(pivot_index[y]), sample->getElement(pairs_index[x][0])) -
                                                       df->operator()(sample->getElement(pivot_index[y]), sample->getElement(pairs_index[x][1])));

                    }

                }

#ifdef ENABLE_DEBUG
                std::cout << "==============================================================================\n" << std::endl;
#endif

                size_t pivotIndex = 0, pos = 0, max_pos = 0;
                double max = std::numeric_limits<double>::min(), max_p = std::numeric_limits<double>::min(), prec = 0.0;

                while (pivotIndex < nPivots)
                {

                    max = std::numeric_limits<double>::lowest();
                    pos = 0;

                    for (size_t x = 0; x < pvt_sz; x++)
                    {

                        if (!bitmap[x])
                        {

                            bitmap[x] = true;
                            prec = 0.0;

                            for (size_t i = 0; i < limInf; i++)
                            {

                                max_p = std::numeric_limits<double>::lowest();
                                max_pos = 0;

                                for (size_t j = 0; j < pvt_sz; j++)
                                {

                                    if((bitmap[j]) && (dist_pivots_pairs[j][i] > max_p))
                                    {

                                        max_p = dist_pivots_pairs[j][i];
                                        max_pos = j;

                                    }

                                }

                                prec += dist_pivots_pairs[max_pos][i]/dist_pairs[i];

                            }

                            prec /= static_cast<double>(limInf);

                            if (prec > max)
                            {

                                max = prec;
                                pos = x;

                            }

                            bitmap[x] = false;

                        }

                    }

                    wrapper->resetIndex(sample->getElement(pivot_index[pos]));

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << pivotIndex << ": " << sample->operator[](pivot_index[pos]) << std::endl;
#endif

                    bitmap[pos] = true;
                    this->setPivot(pivotIndex++, sample->getElement(pivot_index[pos]));

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                bitmap.clear();
                pivot_index.clear();
                pairs_index.clear();
                dist_pairs.clear();
                dist_pivots_pairs.clear();
                wrapper.reset();

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

#endif //GERVLIB_HFIPIVOTS_H
