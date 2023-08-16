//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_WDRPIVOTS_H
#define GERVLIB_WDRPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class WDRPivots : public Pivot<O, T>
    {

    private:
        inline static std::variant<size_t, double> pivotSize = -1.0, limInfSampleSize = (size_t)300;
        inline static size_t numIter = 500, alpha = 2;

    private:
        double func(size_t cand_size, size_t pair_size, size_t p_num, std::vector<double> &dist_pair, std::vector<std::vector<double>> &dist_pair_cand, std::vector<size_t> &pivot_index)
        {

            double sum = 0, _dist_pair, _dist_pair_pivot, max = std::numeric_limits<double>::min();

            for(size_t i = 0; i < pair_size; i++)
            {

                _dist_pair = dist_pair[i];
                max = std::numeric_limits<double>::min();

                for(size_t j = 0; j < p_num; j++)
                {

                    _dist_pair_pivot = dist_pair_cand[i][pivot_index[j]];

                    if(_dist_pair_pivot > max)
                    {

                        max = _dist_pair_pivot;

                    }

                }

                sum += _dist_pair * pow(1-max/_dist_pair, static_cast<double>(alpha))/static_cast<double>(cand_size);

            }

            return sum;

        }

    public:
        WDRPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::WDR); }

        WDRPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                  std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                  size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::WDR);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~WDRPivots() override = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            memcpy(data.get() + offset, &this->numIter, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &this->alpha, sizeof(size_t));
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

            memcpy(&this->numIter, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&this->alpha, _data.get() + offset, sizeof(size_t));
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

            size_t ans = sizeof(size_t) * 2 + //number of iterations and alpha
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
            os << "Number of iterations: " << this->numIter << std::endl;
            os << "Alpha: " << this->alpha << std::endl;
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

        void setNumberOfIterations(size_t n){ numIter = n; }

        size_t getNumberOfIterations(){ return numIter; }

        void setAlpha(size_t a){ alpha = a; }

        size_t getAlpha(){ return alpha; }

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
                    throw std::runtime_error("WDRPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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
                    throw std::invalid_argument("WDRPivots::operator(): Pivot size must be greater than the number of pivots.");

                if(limInfSampleSize.index() == 0)
                    limInf = std::get<size_t>(limInfSampleSize);
                else
                    limInf = std::floor(std::get<double>(limInfSampleSize) * sample->getCardinality());

                if (pvt_sz > sample->getCardinality())
                    pvt_sz = sample->getCardinality();

                std::vector<size_t> cand = utils::generateRandomNumbers(0, sample->getCardinality() - 1, pvt_sz, false, this->getSeed());
                std::vector<std::vector<size_t>> pairs(limInf, std::vector<size_t>(2, 0));
                std::vector<bool> bitmap(pvt_sz, true);
                std::vector<size_t> pivot_index(nPivots, 0);
                std::vector<double> dist_pair(limInf, 0.0);
                std::vector<std::vector<double>> dist_pair_cand(limInf, std::vector<double>(pvt_sz, 0.0));
                size_t p_num = 0, pos = 0, iter = numIter;
                bool flag = true;
                double min = std::numeric_limits<double>::max(), f, tmp;
                utils::Random<size_t> rGen(this->getSeed());

#ifdef ENABLE_DEBUG
                std::cout << "============= Pairs candidates: " << limInf << " ===============================================\n" << std::endl;
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
                std::cout << "=================================================================================\n" << std::endl;

                for(size_t i = 0; i < pvt_sz; i++)
                {

                    std::cout << "Candidate " << i << ": " << sample->operator[](cand[i]) << std::endl;

                }

                std::cout << "=================================================================================\n" << std::endl;
#endif

                for(size_t i = 0; i < limInf; i++)
                {

                    dist_pair[i] = df->operator()(sample->getElement(pairs[i][0]), sample->getElement(pairs[i][1]));

                    for(size_t j = 0; j < pvt_sz; j++)
                    {

                        dist_pair_cand[i][j] = fabs(df->operator()(sample->getElement(pairs[i][0]), sample->getElement(cand[j])) -
                                                    df->operator()(sample->getElement(pairs[i][1]), sample->getElement(cand[j])));

                    }

                }

                while (iter--)
                {

                    flag = true;

                    if (p_num < nPivots)
                    {

                        min = std::numeric_limits<double>::max();

                        for (size_t i = 0; i < pvt_sz; i++)
                        {

                            if (bitmap[i])
                            {

                                pivot_index[p_num] = i;
                                f = func(pvt_sz, limInf, p_num + 1, dist_pair, dist_pair_cand, pivot_index);

                                if (f < min)
                                {

                                    min = f;
                                    pos = i;

                                }

                            }

                        }

                        pivot_index[p_num] = pos;
                        bitmap[pos] = false;
                        p_num++;

                    }

                    tmp = func(pvt_sz, limInf, p_num, dist_pair, dist_pair_cand, pivot_index);

                    for (size_t i = 0; i < p_num; i++)
                    {

                        pos = pivot_index[i];

                        for (size_t j = 0; j < pvt_sz; j++)
                        {

                            if (bitmap[j])
                            {

                                pivot_index[i] = j;
                                f = func(pvt_sz, limInf, p_num, dist_pair, dist_pair_cand, pivot_index);

                                if (tmp > f)
                                {

                                    bitmap[j] = false;
                                    tmp = f;
                                    pos = j;
                                    flag = false;

                                }

                            }

                        }

                        pivot_index[i] = pos;

                    }

                    if (flag && p_num == nPivots)
                        break;

                }

                for (size_t i = 0; i < p_num; i++)
                {

                    this->setPivot(i, sample->getElement(cand[pivot_index[i]]));

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << i << ": " << sample->operator[](cand[pivot_index[i]]) << std::endl;
#endif

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                cand.clear();
                pairs.clear();
                bitmap.clear();
                pivot_index.clear();
                dist_pair.clear();
                dist_pair_cand.clear();

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

#endif //GERVLIB_WDRPIVOTS_H
