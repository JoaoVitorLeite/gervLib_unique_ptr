//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_SELECTIONPIVOTS_H
#define GERVLIB_SELECTIONPIVOTS_H

#include "Pivot.h"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class SelectionPivots : public Pivot<O, T>
    {

    private:
        inline static std::variant<size_t, double> nGroups = 0.5;

    public:
        SelectionPivots(): Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::SELECTION); }

        SelectionPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                        size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::SELECTION);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~SelectionPivots() override = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            size = this->nGroups.index();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (this->nGroups.index() == 0)
            {

                memcpy(data.get() + offset, &std::get<size_t>(this->nGroups), sizeof(size_t));
                offset += sizeof(size_t);

            }
            else
            {

                memcpy(data.get() + offset, &std::get<double>(this->nGroups), sizeof(double));
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
                this->nGroups = sz2;

            }
            else
            {

                double sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                this->nGroups = sz2;

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

            size_t ans = sizeof(size_t) + (this->nGroups.index() == 0 ? sizeof(size_t) : sizeof(double)) + //number of groups
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
            os << "Number of groups: " << configure::variant2string(this->nGroups) << std::endl;
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Number of drop pivots: " << this->drop << std::endl;
            os << "Pivot Number of Pivots: " << this->pivots->getCardinality() << std::endl;
            os << "Pivot Pivots: " << std::endl;
            os << *this->pivots << std::endl;

        }

        std::variant<size_t, double> getNumberOfGroups(){ return nGroups; }

        void setNumberOfGroups(size_t size){ nGroups = size; }

        void setNumberOfGroups(double size){ nGroups = size; }

        void setNumberOfGroups(std::variant<size_t, double> size){ nGroups = size; }

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
                    throw std::runtime_error("SelectionPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
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

                size_t sz = 0;

                if(nGroups.index() == 0)
                    sz = std::get<size_t>(nGroups);
                else
                    sz = std::floor(std::get<double>(nGroups) * sample->getCardinality());

                std::vector<size_t> pvt_group(nPivots);
                double max = std::numeric_limits<double>::min(), len = static_cast<double>((nPivots * (nPivots - 1)))/2.0, dist, mean;
                utils::Random rGen(this->getSeed());

                for (size_t g = 0; g < sz; g++)
                {

                    std::vector<size_t> indexes = utils::generateRandomNumbers(0, sample->getCardinality()-1, nPivots, false,
                                                                               rGen());

                    mean = 0.0;

                    for (size_t i = 0; i < nPivots; i++) {

                        for (size_t j = i + 1; j < nPivots; j++){

                            dist = df->operator()(sample->getElement(indexes[i]), sample->getElement(indexes[j]));
                            mean += dist/len;

                        }

                    }

                    if (mean > max)
                    {

                        max = mean;
                        std::copy(indexes.begin(), indexes.end(), pvt_group.begin());

                    }

#ifdef ENABLE_DEBUG
                    std::cout << "Group " << g << ": " << mean << std::endl;

                    for (size_t z = 0; z < nPivots; z++)
                    {

                        std::cout << "Element: " << sample->operator[](indexes[z]) << std::endl;

                    }
#endif

                }

                for (size_t i = 0; i < nPivots; i++) {

                    this->setPivot(i, sample->getElement(pvt_group[i]));

                }


                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                pvt_group.clear();

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

#endif //GERVLIB_SELECTIONPIVOTS_H
