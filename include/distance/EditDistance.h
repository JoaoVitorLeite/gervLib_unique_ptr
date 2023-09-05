//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_EDITDISTANCE_H
#define GERVLIB_EDITDISTANCE_H

#include "DistanceFunction.h"

namespace gervLib::distance
{
    template <typename T>
    requires HasSize<T> && IsChar<T>
    class EditDistance: public DistanceFunction<T> {

    public:
        explicit EditDistance(): DistanceFunction<T>() { this->setDistanceType(LEVENSHTEIN); }
        ~EditDistance() override = default;

        double operator()(const T& a, const T& b) override {

            if(a.size() != b.size())
                throw std::invalid_argument("Edit distance requires objects of the same size.");

            this->incrementStatistics();
            size_t s = a.size();
            double dist = 0.0;

            for (size_t x = 0; x < s; x++) {

                size_t s1 = a[x].size();
                size_t s2 = b[x].size();

                auto* column = new size_t[s1 + 1];
                std::iota(column, column + s1 + 1, 0);

                for(size_t i = 1; i <= s2; i++) {
                    column[0] = i;
                    size_t lastDiagonal = i - 1;
                    for(size_t j = 1; j <= s1; j++) {
                        size_t oldDiagonal = column[j];
                        column[j] = std::min({column[j] + 1, column[j - 1] + 1, lastDiagonal + (a[x][j - 1] == b[x][i - 1] ? 0 : 1)});
                        lastDiagonal = oldDiagonal;
                    }
                }

                dist += static_cast<double>(column[s1]);
                delete[] column;

            }

            return dist;

        }

    };
}

#endif //GERVLIB_EDITDISTANCE_H
