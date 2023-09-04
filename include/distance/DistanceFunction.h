//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_DISTANCEFUNCTION_H
#define GERVLIB_DISTANCEFUNCTION_H

#include <concepts>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <iostream>
#include <vector>
#include <cstring>
#include "Serialize.h"

namespace gervLib::distance {

    enum DISTANCE_TYPE {
        EUCLIDEAN,
        LEVENSHTEIN,
        UNKNOWN
    };

    std::map<DISTANCE_TYPE, std::string> distanceTypeMap = {
            {EUCLIDEAN,   "EUCLIDEAN"},
            {LEVENSHTEIN, "LEVENSHTEIN"}
    };

    std::map<std::string, DISTANCE_TYPE> distanceNameMap = {
            {"EUCLIDEAN",   EUCLIDEAN},
            {"LEVENSHTEIN", LEVENSHTEIN}
    };

    template<typename T>
    concept IsDouble = requires(T t) {
        { t[0] } -> std::convertible_to<double>;
    };

    template<typename T>
    concept IsChar = requires(T t) {
        { t[0] } -> std::convertible_to<std::vector<char>>;
    };

    template<typename T>
    concept HasSize = requires(T t) {
        { t.size() } -> std::convertible_to<size_t>;
    };

    template<typename T> requires HasSize<T> && (IsDouble<T> || IsChar<T>)
    class DistanceFunction : public gervLib::serialize::Serialize {

    private:
        inline static size_t distCount = 0;

    protected:
        DISTANCE_TYPE distanceType;

    protected:
        void incrementStatistics() {
            distCount++;
        }

    public:
        explicit DistanceFunction() : distanceType(UNKNOWN) { }

        virtual ~DistanceFunction() = default;

        [[nodiscard]] size_t getDistanceCount() const {
            return distCount;
        }

        void resetStatistics() {
            distCount = 0;
        }

        DistanceFunction<T> &operator=(const DistanceFunction<T> &other) {
            distCount = other.distCount;
            return *this;
        }

        virtual double operator()(const T &a, const T &b) = 0;

        double getDistance(const T &a, const T &b) {
            return operator()(a, b);
        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0;

            memcpy(data.get() + offset, &distCount, sizeof(size_t));
            offset += sizeof(size_t);

            return data;
        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {
            size_t offset = 0;

            memcpy(&distCount, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            _data.reset();

        }

        size_t getSerializedSize() override {

            return sizeof(size_t);

        }

        [[nodiscard]] std::string getDistanceFunctionName() const { return distanceTypeMap[distanceType]; }

        bool isEqual(std::unique_ptr<DistanceFunction<T>>& other) {
            return ((distCount == other->distCount) && (getDistanceFunctionName() == other->getDistanceFunctionName()));
        }

        bool operator==(DistanceFunction<T> *other) {
            return isEqual(other);
        }

        bool operator!=(DistanceFunction<T> *other) {
            return !isEqual(other);
        }

        void setDistanceType(DISTANCE_TYPE _distanceType) {
            this->distanceType = _distanceType;
        }

        [[nodiscard]] DISTANCE_TYPE getDistanceType() const {
            return distanceType;
        }

        friend std::ostream &operator<<(std::ostream &os, const DistanceFunction<T> &distanceFunction) {
            os << "DistanceFunction: " << distanceFunction.getDistanceFunctionName() << std::endl;
            os << "Distance count: " << distanceFunction.getDistanceCount() << std::endl;
            return os;
        }


    };

}


#endif //GERVLIB_DISTANCEFUNCTION_H
