//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_EUCLIDEANDISTANCE_H
#define GERVLIB_EUCLIDEANDISTANCE_H

#include "DistanceFunction.h"

namespace gervLib::distance
{
    template <typename T>
    requires HasSize<T> && IsDouble<T>
    class EuclideanDistance: public DistanceFunction<T> {

    public:
        explicit EuclideanDistance(size_t distance = 0): DistanceFunction<T>(distance) { this->setDistanceType(EUCLIDEAN); }
        ~EuclideanDistance() override = default;

        double operator()(const T& a, const T& b) override {

            if(a.size() != b.size())
                throw std::invalid_argument("Euclidean distance requires objects of the same size.");

            this->incrementStatistics();
            double sum = 0;
            for (size_t i = 0; i < a.size(); i++) {
                sum += (a[i] - b[i]) * (a[i] - b[i]);
            }
            return sqrt(sum);
        }

    };

}

#endif //GERVLIB_EUCLIDEANDISTANCE_H
