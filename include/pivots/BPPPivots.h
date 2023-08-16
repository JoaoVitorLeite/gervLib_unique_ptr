//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_BPPPIVOTS_H
#define GERVLIB_BPPPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class BPPPivots : public Pivot<O, T>
    {

    private:
        inline static std::variant<size_t, double> pivotSize = -1.0, limInfSampleSize = (size_t)300;

    public:
        BPPPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::BPP); }

        BPPPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                  std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                  size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::BPP);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~BPPPivots() override = default;

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
                    throw std::runtime_error("BPPPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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
                    throw std::invalid_argument("BPPPivots::operator(): Pivot size must be greater than the number of pivots.");

                if(limInfSampleSize.index() == 0)
                    limInf = std::get<size_t>(limInfSampleSize);
                else
                    limInf = std::floor(std::get<double>(limInfSampleSize) * sample->getCardinality());

                if (pvt_sz > sample->getCardinality())
                    pvt_sz = sample->getCardinality();

                if (limInf > sample->getCardinality())
                    limInf = sample->getCardinality();

                size_t cand = pvt_sz, pos = 0;
                double d = 0.0, min = std::numeric_limits<double>::max();
                utils::Random<size_t> rGen(this->getSeed());
                std::vector<size_t> pivot_index = utils::generateRandomNumbers(0, sample->getCardinality() - 1, pvt_sz, false, rGen());
                std::vector<size_t> aux_index = utils::generateRandomNumbers(0, sample->getCardinality() - 1, limInf, false, rGen());
                std::vector<bool> bitmap(pvt_sz, true);
                std::vector<std::vector<double>> dist(pvt_sz, std::vector<double>(limInf, 0.0));
                std::vector<std::vector<double>> pr_mat(pvt_sz, std::vector<double>(pvt_sz+1, 0.0));
                std::vector<std::vector<size_t>> rank(pvt_sz, std::vector<size_t>(limInf, 1));

#ifdef ENABLE_DEBUG

                std::cout << "==============================================================================================" << std::endl;

                for (size_t i = 0; i < pvt_sz; i++)
                    std::cout << "Pivot candidates " << i << ": " << sample->operator[](pivot_index[i]) << std::endl;

                std::cout << std::endl << std::endl;

                for (size_t i = 0; i < limInf; i++)
                    std::cout << "Auxiliary " << i << ": " << sample->operator[](aux_index[i]) << std::endl;

                std::cout << "==============================================================================================\n" << std::endl;

#endif

                for(size_t i = 0; i < pvt_sz; i++)
                    for(size_t j = 0; j < limInf; j++)
                        dist[i][j] = df->operator()(sample->getElement(aux_index[j]), sample->getElement(pivot_index[i]));

                for(size_t i = 0; i < pvt_sz; i++)
                {

                    for(size_t j = 0; j < limInf; j++)
                    {

                        d = dist[i][j];

                        for(size_t k = 0; k < pvt_sz; k++)
                        {

                            if(k != i && dist[k][j] < d)
                                rank[i][j]++;

                        }

                    }

                }

                while(cand > nPivots)
                {

                    min = std::numeric_limits<double>::max();
                    pos = 0;

                    for(size_t i = 0; i < pvt_sz; i++)
                    {

                        if(bitmap[i])
                        {

                            bitmap[i] = false;

                            double std = 0.0;

                            for(size_t j = 0; j < pvt_sz; j++)
                            {

                                if(bitmap[j])
                                {

                                    size_t k = 1, pr = 0, sum = 0;

                                    while(k < cand)
                                    {

                                        for(size_t z = 0; z < limInf; z++)
                                        {

                                            if(rank[j][z] == k) pr++;

                                        }

                                        sum += pr;
                                        pr_mat[j][k-1] = static_cast<double>(pr);
                                        pr = 0;
                                        k++;

                                    }

                                    pr_mat[j][cand] = (static_cast<double>(sum)*1.0)/(static_cast<double>(cand)-1);

                                    for(size_t x = 0; x < cand-1; x++)
                                    {

                                        std += (pr_mat[j][x] - pr_mat[j][cand])*(pr_mat[j][x] - pr_mat[j][cand]);

                                    }

                                }

                            }

                            std /= (static_cast<double>(cand)-1)*(static_cast<double>(cand)-1);

                            if(std < min)
                            {

                                min = std;
                                pos = i;

                            }

                            bitmap[i] = true;

                        }

                    }

                    bitmap[pos] = false;
                    cand--;

                }

                pos = 0;

                for(size_t i = 0; i < pvt_sz; i++)
                {

                    if(bitmap[i])
                    {

#ifdef ENABLE_DEBUG
                        std::cout << "Pivot " << pos << ": " << sample->operator[](pivot_index[i]) << std::endl;
#endif
                        this->setPivot(pos++, sample->getElement(pivot_index[i]));

                    }

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                pivot_index.clear();
                aux_index.clear();
                bitmap.clear();
                rank.clear();
                pr_mat.clear();
                dist.clear();

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

#endif //GERVLIB_BPPPIVOTS_H
