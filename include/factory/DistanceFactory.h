//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_DISTANCEFACTORY_H
#define GERVLIB_DISTANCEFACTORY_H

#include "EuclideanDistance.h"
#include "EditDistance.h"
#include <memory>

namespace gervLib::distance {

    template <typename T>
    requires HasSize<T> && (IsDouble<T> || IsChar<T>)
    class DistanceFactory {

    public:
        static std::unique_ptr<DistanceFunction<T>> createDistanceFunction(DISTANCE_TYPE distanceType) {

            if (distanceType == EUCLIDEAN)
            {
                if constexpr (IsDouble<T>)
                    return std::make_unique<EuclideanDistance<T>>();
                else
                    throw std::invalid_argument("EditDistance requires IsChar<T>");
            }
            else if(distanceType == LEVENSHTEIN)
            {
                if constexpr (IsChar<T>)
                    return std::make_unique<EditDistance<T>>();
                else
                    throw std::invalid_argument("EditDistance requires IsChar<T>");
            }
            else
                throw std::invalid_argument("Invalid DISTANCE_TYPE");

        }


    };

}

#endif //GERVLIB_DISTANCEFACTORY_H
