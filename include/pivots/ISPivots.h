//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_ISPIVOTS_H
#define GERVLIB_ISPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class ISPivots : public Pivot<O, T>
    {

    private:
        inline static std::variant<size_t, double> pivotSize = -1.0, limInfSampleSize = (size_t)300;

    public:
        ISPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::IS); }

        ISPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                 std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                 size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::IS);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~ISPivots() override = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

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

            size_t ans = sizeof(size_t) + (this->pivotSize.index() == 0 ? sizeof(size_t) : sizeof(double)) + //pivot size
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
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Number of drop pivots: " << this->drop << std::endl;
            os << "Pivot Number of Pivots: " << this->pivots->getCardinality() << std::endl;
            os << "Pivot Pivots: " << std::endl;
            os << *this->pivots << std::endl;

        }

        void setPivotSize(size_t sz){ pivotSize = sz; }

        void setPivotSize(double sz){ pivotSize = sz; }

        std::variant<size_t, double> getPivotSize(){ return pivotSize; }

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
                    throw std::runtime_error("ISPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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
                    throw std::invalid_argument("ISPivots::operator(): Pivot size must be greater than the number of pivots.");

                if(limInfSampleSize.index() == 0)
                    limInf = std::get<size_t>(limInfSampleSize);
                else
                    limInf = std::floor(std::get<double>(limInfSampleSize) * sample->getCardinality());

                if (pvt_sz > sample->getCardinality())
                    pvt_sz = sample->getCardinality();

                std::vector<bool> bitmap(pvt_sz, true);
                std::vector<size_t> candidates = utils::generateRandomNumbers(0, sample->getCardinality() - 1, pvt_sz, false, this->getSeed());
                std::vector<std::vector<size_t>> pairs(limInf, std::vector<size_t>(2, 0));
                std::vector<double> svg(limInf, 0.0);
                double max = std::numeric_limits<double>::min(), q = 0.0, dfs = 0.0, d = 0.0;
                size_t pos = 0;
                utils::Random<size_t> rGen(this->getSeed());

#ifdef ENABLE_DEBUG
                std::cout << "==============================================================================\n" << std::endl;

                for(size_t i = 0; i < pvt_sz; i++)
                    std::cout << "Candidate " << i << ": " << sample->operator[](candidates[i]) << std::endl;

                std::cout << "============= Pairs: " << limInf << " =======================================================\n" << std::endl;
#endif

                for(size_t x = 0; x < limInf; x++)
                {

                    while (pairs[x][0] == pairs[x][1])
                    {

                        pairs[x][0] = rGen(0, sample->getCardinality()-1);
                        pairs[x][1] = rGen(0, sample->getCardinality()-1);

                    }

#ifdef ENABLE_DEBUG
                    std::cout << "Pair " << x << ": " << sample->operator[](pairs[x][0]) << " <->  " << sample->operator[](pairs[x][1]) << std::endl;
#endif

                }

#ifdef ENABLE_DEBUG
                std::cout << "==============================================================================\n" << std::endl;
#endif

                for (size_t x = 0; x < nPivots; x++)
                {

                    max = std::numeric_limits<double>::min();
                    pos = 0;

                    for (size_t y = 0; y < pvt_sz; y++)
                    {

                        if (bitmap[y])
                        {

                            dfs = 0.0;
                            d = 0.0;

                            for (size_t z = 0; z < limInf; z++)
                            {

                                d = fabs(df->operator()(sample->getElement(candidates[y]), sample->getElement(pairs[z][1])) -
                                         df->operator()(sample->getElement(candidates[y]), sample->getElement(pairs[z][0])));

                                if (svg[z] > d)
                                    dfs += svg[z];
                                else
                                    dfs += d;

                            }

                            if (dfs > max)
                            {

                                max = dfs;
                                pos = y;

                            }

                        }


                    }

                    bitmap[pos] = false;

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << x << ": " << sample->operator[](candidates[pos]) << std::endl;
#endif

                    this->setPivot(x, sample->getElement(candidates[pos]));

                    for (size_t i = 0; i < limInf; i++)
                    {

                        q = fabs(df->operator()(sample->getElement(candidates[pos]), sample->getElement(pairs[i][1])) -
                                 df->operator()(sample->getElement(candidates[pos]), sample->getElement(pairs[i][0])));

                        if (q > svg[i])
                            svg[i] = q;

                    }


                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                bitmap.clear();
                candidates.clear();
                pairs.clear();
                svg.clear();

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

#endif //GERVLIB_ISPIVOTS_H
