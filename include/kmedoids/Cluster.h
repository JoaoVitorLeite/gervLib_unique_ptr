//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_CLUSTER_H
#define GERVLIB_CLUSTER_H

#include <vector>
#include "Dataset.h"
#include "DistanceFunction.h"

namespace gervLib::kmedoids {

    template<typename O, typename T>
    class Cluster {

    private:
        std::vector<dataset::BasicArrayObject<O, T>> objects;
        dataset::BasicArrayObject<O, T> medoid;

    public:
        Cluster() = default;

        explicit Cluster(dataset::BasicArrayObject<O, T> _medoid) : medoid(_medoid) {}

        ~Cluster() { clear(); }

        void add(dataset::BasicArrayObject<O, T> obj) { objects.push_back(obj); }

        dataset::BasicArrayObject<O, T>& getMedoid() { return medoid; }

        dataset::BasicArrayObject<O, T>& get(size_t pos) {

            utils::check_range(0, objects.size(), pos, "Cluster::get");
            return objects[pos];

        }

        void setMedoid(dataset::BasicArrayObject<O, T> _medoid) { medoid = _medoid; }

        void clear() { objects.clear(); }

        size_t size() { return objects.size(); }

        friend std::ostream &operator<<(std::ostream &os, const Cluster &cluster) {

            os << "Medoid: " << cluster.medoid << std::endl;
            os << "Objects: " << std::endl;

            for (auto &object: cluster.objects) {
                os << object << std::endl;
            }

            return os;

        }

    };

}

#endif //GERVLIB_CLUSTER_H
