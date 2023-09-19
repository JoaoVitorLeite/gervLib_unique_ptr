//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_DISTANCEFACTORY_H
#define GERVLIB_DISTANCEFACTORY_H

#include "EuclideanDistance.h"
//#include "EditDistance.h"
#include <memory>

namespace gervLib::distance {

    template <typename T>
    class DistanceFactory {

    public:
        static std::unique_ptr<DistanceFunction<T>> createDistanceFunction(DISTANCE_TYPE distanceType) {

            if (distanceType == EUCLIDEAN)
            {
                return std::make_unique<EuclideanDistance<T>>();
            }
            else
                throw std::invalid_argument("Invalid DISTANCE_TYPE");

        }


    };

}

#endif //GERVLIB_DISTANCEFACTORY_H
