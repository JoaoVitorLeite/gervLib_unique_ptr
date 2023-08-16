//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_PIVOT_H
#define GERVLIB_PIVOT_H

#include "DistanceFunction.h"
#include "Dataset.h"

namespace gervLib::pivots {

    enum PIVOT_TYPE {
        RANDOM, CONVEX, KMEDOIDS, MAXSEPARATED, MAXVARIANCE, SELECTION, PCA, SSS, FFT, HFI, IS, WDR, BPP
    };

    std::map<PIVOT_TYPE, std::string> PIVOT_TYPE2STR = {
            {PIVOT_TYPE::RANDOM,       "RANDOM"},
            {PIVOT_TYPE::CONVEX,       "CONVEX"},
            {PIVOT_TYPE::KMEDOIDS,     "KMEDOIDS"},
            {PIVOT_TYPE::MAXSEPARATED, "MAXSEPARATED"},
            {PIVOT_TYPE::MAXVARIANCE,  "MAXVARIANCE"},
            {PIVOT_TYPE::SELECTION,    "SELECTION"},
            {PIVOT_TYPE::PCA,          "PCA"},
            {PIVOT_TYPE::SSS,          "SSS"},
            {PIVOT_TYPE::FFT,          "FFT"},
            {PIVOT_TYPE::HFI,          "HFI"},
            {PIVOT_TYPE::IS,           "IS"},
            {PIVOT_TYPE::WDR,          "WDR"},
            {PIVOT_TYPE::BPP,          "BPP"}
    };

    std::map<std::string, PIVOT_TYPE> STR2PIVOT_TYPE = {
            {"RANDOM",       PIVOT_TYPE::RANDOM},
            {"CONVEX",       PIVOT_TYPE::CONVEX},
            {"KMEDOIDS",     PIVOT_TYPE::KMEDOIDS},
            {"MAXSEPARATED", PIVOT_TYPE::MAXSEPARATED},
            {"MAXVARIANCE",  PIVOT_TYPE::MAXVARIANCE},
            {"SELECTION",    PIVOT_TYPE::SELECTION},
            {"PCA",          PIVOT_TYPE::PCA},
            {"SSS",          PIVOT_TYPE::SSS},
            {"FFT",          PIVOT_TYPE::FFT},
            {"HFI",          PIVOT_TYPE::HFI},
            {"IS",           PIVOT_TYPE::IS},
            {"WDR",          PIVOT_TYPE::WDR},
            {"BPP",          PIVOT_TYPE::BPP}
    };

    template <typename O, typename T>
    class Pivot : public serialize::Serialize
    {

    protected:
        std::unique_ptr<dataset::Dataset<O, T>> pivots;
        size_t seed;
        PIVOT_TYPE type;
        long long elapsedTime;
        inline static std::variant<size_t, double> sampleSize = 1.0;
        inline static size_t drop = 0;

    protected:
        void incrementElapsedTime(long long time)
        {

            elapsedTime += time;

        }

        bool needSample(std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {

            if (sampleSize.index() == 0)
            {

                return dataset->getCardinality() != std::get<size_t>(sampleSize);

            }
            else
            {

                return std::get<double>(sampleSize) != 1.0;

            }

        }

    public:
        // Constructor and destructor

        Pivot() : pivots(std::make_unique<dataset::Dataset<O, T>>()), seed(0), type(PIVOT_TYPE::RANDOM), elapsedTime(0) {}

        virtual ~Pivot()
        {
            if (pivots != nullptr) {
                pivots->clear();
                pivots.reset();
            }
        }

        // Virtual methods

        virtual void print(std::ostream& os) const
        {

            os << "Pivot Type: " << type << std::endl;
            os << "Pivot Seed: " << seed << std::endl;
            os << "Pivot Elapsed Time: " << elapsedTime << std::endl;
            os << "Pivot Sample Size: " << configure::variant2string(this->sampleSize) << std::endl;
            os << "Number of drop pivots: " << drop << std::endl;

            if (pivots == nullptr)
                os << "Pivot Number of Pivots: 0" << std::endl;
            else
                os << "Pivot Number of Pivots: " << pivots->getCardinality() << std::endl;

            os << "Pivot Pivots: " << std::endl;
            if (pivots != nullptr)
                os << *pivots << std::endl;
            else
                os << "Pivot Pivots: NULL" << std::endl;

        }

        virtual void operator()(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                                std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots) { }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            memcpy(data.get() + offset, &seed, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &drop, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &elapsedTime, sizeof(long long));
            offset += sizeof(long long);

            size = sampleSize.index();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (sampleSize.index() == 0)
            {

                memcpy(data.get() + offset, &std::get<size_t>(sampleSize), sizeof(size_t));
                offset += sizeof(size_t);

            }
            else
            {

                memcpy(data.get() + offset, &std::get<double>(sampleSize), sizeof(double));
                offset += sizeof(double);

            }

            size = pivots == nullptr ? 0 : pivots->getSerializedSize();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            if (pivots != nullptr)
            {

                std::unique_ptr<u_char[]> aux = pivots->serialize();
                memcpy(data.get() + offset, aux.get(), size);
                offset += size;
                aux.reset();

            }

            size = PIVOT_TYPE2STR[type].size();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, PIVOT_TYPE2STR[type].c_str(), size);
            offset += size;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, size;

            memcpy(&seed, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&drop, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&elapsedTime, _data.get() + offset, sizeof(long long));
            offset += sizeof(long long);

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                size_t sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);
                sampleSize = sz2;

            }
            else
            {

                double sz2;
                memcpy(&sz2, _data.get() + offset, sizeof(double));
                offset += sizeof(double);
                sampleSize = sz2;

            }

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {

                clear();

            }
            else
            {

                pivots->clear();
                pivots = std::make_unique<dataset::Dataset<O, T>>();
                std::unique_ptr<u_char[]> aux(new u_char[size]);
                memcpy(aux.get(), _data.get() + offset, size);
                offset += size;
                pivots->deserialize(std::move(aux));

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

            size_t ans = sizeof(size_t) + //seed
                         sizeof(size_t) + //drop
                         sizeof(long long) + //time
                         sizeof(size_t) + (sampleSize.index() == 0 ? sizeof(size_t) : sizeof(double)); //number to indicate variant possibility and variable size

            ans += sizeof(size_t) + (pivots == nullptr ? 0 : pivots->getSerializedSize());
            std::string aux = PIVOT_TYPE2STR[getPivotType()];
            ans += sizeof(size_t) + aux.size();

            return ans;

        }

        // Setters

        void setSeed(size_t _seed)
        {

            seed = _seed;

        }

        void setPivot(size_t pos, dataset::BasicArrayObject<O, T>& value)
        {

            utils::check_range(0, getNumberOfPivots()-1, pos, "Pivot::setPivot");
            pivots->operator[](pos) = value;

        }

        void setPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {

            clear();
            pivots = std::make_unique<dataset::Dataset<O, T>>(*dataset);

        }

        void setPivotType(PIVOT_TYPE _type)
        {

            type = _type;

        }

        void setNumberOfDropPivots(size_t _drop)
        {

            drop = _drop;

        }

        void setSampleSize(size_t sz)
        {

            sampleSize = sz;

        }

        void setSampleSize(double sz)
        {

            sampleSize = sz;

        }

        void setSampleSize(std::variant<size_t, double> sz)
        {

            sampleSize = sz;

        }

        void setElapsedTime(long long time)
        {

            elapsedTime = time;

        }

        void setNumberOfPivots(size_t nPivots)
        {

            if (pivots == nullptr) {
                pivots = std::make_unique<dataset::Dataset<O, T>>();
                pivots->setData(std::vector<dataset::BasicArrayObject<O, T>>(nPivots));
            }
            else {
                pivots->setData(std::vector<dataset::BasicArrayObject<O, T>>());
                pivots->setData(std::vector<dataset::BasicArrayObject<O, T>>(nPivots));
            }
        }

        // Getters

        dataset::BasicArrayObject<O, T>& getPivot(size_t pos)
        {

            if (pivots == nullptr)
                throw std::runtime_error("Pivot::getPivot: pivots is null");

            utils::check_range(0, pivots->getCardinality()-1, pos, "Pivot::getPivot");
            return pivots->operator[](pos);

        }

        std::unique_ptr<dataset::Dataset<O, T>>& getPivots()
        {

            return pivots;

        }

        PIVOT_TYPE getPivotType()
        {

            return type;

        }

        size_t getNumberOfDropPivots()
        {

            return drop;

        }

        size_t getNumberOfPivots()
        {

            if (pivots == nullptr)
                return 0;
            else
                return pivots->getCardinality();

        }

        std::variant<size_t, double> getSampleSize()
        {

            return sampleSize;

        }

        size_t getSeed()
        {

            return seed;

        }

        long long getElapsedTime()
        {

            return elapsedTime;

        }

        // Operators

        bool operator==(std::unique_ptr<Pivot<O, T>> other)
        {

            return isEqual(other);

        }

        bool operator!=(std::unique_ptr<Pivot<O, T>> other)
        {

            return !isEqual(other);

        }

        // Methods

        void clear()
        {
            if (pivots != nullptr) {
                pivots->clear();
                pivots.reset();
            }
        }

        bool isEqual(std::unique_ptr<Pivot<O, T>>& other)
        {

            if (pivots != nullptr && other->pivots == nullptr)
                return false;

            if (pivots == nullptr && other->pivots != nullptr)
                return false;

            if (pivots != nullptr && other->pivots != nullptr)
                if (pivots->getCardinality() != other->pivots->getCardinality())
                    return false;

            if(type != other->type || seed != other->seed)
                return false;

            for(size_t i = 0; i < getNumberOfPivots(); i++)
                if(pivots->operator[](i) != other->getPivot(i))
                    return false;

            return true;

        }

        void generatePivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots)
        {

            operator()(dataset, df, nPivots);

        }

    };

    template <typename O, typename T>
    std::ostream& operator<<(std::ostream& os, const Pivot<O, T>& printable) {
        printable.print(os);
        return os;
    }

}

#endif //GERVLIB_PIVOT_H
