//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_KMEDOIDSPIVOTS_H
#define GERVLIB_KMEDOIDSPIVOTS_H

#include "Pivot.h"
#include "Kmedoids.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class KmedoidsPivots : public Pivot<O, T>
    {

    private:
        inline static size_t numberOfIterations = 500;

    public:
        KmedoidsPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::KMEDOIDS); }

        KmedoidsPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                       std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                       size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::FFT);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~KmedoidsPivots() override = default;

        void print(std::ostream& os) const override
        {

            os << "Pivot Type: " << this->type << std::endl;
            os << "Pivot Seed: " << this->seed << std::endl;
            os << "Pivot Elapsed Time: " << this->elapsedTime << std::endl;
            os << "Number of Iterations: " << this->numberOfIterations << std::endl;
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Number of drop pivots: " << this->drop << std::endl;
            os << "Pivot Number of Pivots: " << this->pivots->getCardinality() << std::endl;
            os << "Pivot Pivots: " << std::endl;
            os << *this->pivots << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            memcpy(data.get() + offset, &this->numberOfIterations, sizeof(size_t));
            offset += sizeof(size_t);

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

            memcpy(&this->numberOfIterations, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

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

            size_t ans = sizeof(size_t) + //number of iterations
                         sizeof(size_t) + //seed
                         sizeof(size_t) + //drop
                         sizeof(long long) + //time
                         sizeof(size_t) + (this->sampleSize.index() == 0 ? sizeof(size_t) : sizeof(double)); //number to indicate variant possibility and variable size

            ans += sizeof(size_t) + (this->pivots == nullptr ? 0 : this->pivots->getSerializedSize());
            std::string aux = PIVOT_TYPE2STR[this->getPivotType()];
            ans += sizeof(size_t) + aux.size();

            return ans;

        }

        void setNumberOfIterations(size_t _numberOfIterations) { numberOfIterations = _numberOfIterations; }

        size_t getNumberOfIterations() { return KmedoidsPivots::numberOfIterations; }

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

                kmedoids::Kmedoids<O, T> aux = kmedoids::Kmedoids<O, T>();
                aux.run(sample, df, nPivots, numberOfIterations, this->getSeed());

                std::vector<dataset::BasicArrayObject<O, T>> pivots = aux.getMedoids();

                for (size_t i = 0; i < nPivots; i++)
                {

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << i << ": " << pivots[i] << std::endl;
#endif

                    this->setPivot(i, pivots[i]);

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

#endif //GERVLIB_KMEDOIDSPIVOTS_H
