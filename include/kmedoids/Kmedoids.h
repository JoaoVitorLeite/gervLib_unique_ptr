//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_KMEDOIDS_H
#define GERVLIB_KMEDOIDS_H

#include "Cluster.h"

namespace gervLib::kmedoids
{
    template<typename O, typename T>
    class Kmedoids {

    private:
        std::vector<Cluster<O, T>> clusters;

        void initialize(std::unique_ptr<dataset::Dataset<O, T>>& dataset, size_t nClusters, size_t seed)
        {

            utils::Random<size_t> random(seed);
            std::set<size_t> medoids;

            while (medoids.size() < nClusters) {
                medoids.insert(random(0, dataset->getCardinality() - 1));
            }

            for (auto &medoid : medoids)
            {

#ifdef ENABLE_DEBUG
                std::cout << "Random Medoid: " << dataset->operator[](medoid) << std::endl;
#endif

                clusters.emplace_back(dataset->getElement(medoid));

            }

        }

        void assignment(std::unique_ptr<dataset::Dataset<O, T>> &dataset,
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> &df)
        {

            for (size_t i = 0; i < clusters.size(); i++)
                clusters[i].clear();

            size_t pos = 0;
            double min = std::numeric_limits<double>::max(), dist;

            for (size_t i = 0; i < dataset->getCardinality(); i++)
            {

                min = std::numeric_limits<double>::max();
                pos = 0;

                for (size_t j = 0; j < clusters.size(); j++)
                {

                    dist = df->operator()(dataset->getElement(i), clusters[j].getMedoid());

                    if (dist < min) {
                        min = dist;
                        pos = j;
                    }

                }

                clusters[pos].add(dataset->getElement(i));

            }

        }

        void reCenter(std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df)
        {

            double min = std::numeric_limits<double>::max(), dist = 0.0, total = 0.0;
            size_t pos = 0, sz = 0;

            for (size_t i = 0; i < clusters.size(); i++)
            {

                min = std::numeric_limits<double>::max(), total = 0;
                pos = 0, sz = clusters[i].size();

                for (size_t j = 0; j < sz; j++)
                {

                    total = 0.0;

                    for (size_t k = 0; k < sz; k++)
                    {

                        dist = df->operator()(clusters[i].get(j), clusters[i].get(k));
                        total += dist/(static_cast<double>(sz)-1.0);

                    }

                    if (total < min) {
                        min = total;
                        pos = j;
                    }

                }

                clusters[i].setMedoid(clusters[i].get(pos));

            }

        }

    public:
        Kmedoids() = default;

        ~Kmedoids()
        {

            for (auto &cluster : clusters) {
                cluster.clear();
            }

        }

        size_t getNumberOfClusters() { return clusters.size(); }

        std::vector<O> getMedoidsOID()
        {

            std::vector<O> medoidsOIDs;

            for (size_t i = 0; i < clusters.size(); i++)
                medoidsOIDs.push_back(clusters[i].getMedoid().getOID());


            return medoidsOIDs;

        }

        std::vector<dataset::BasicArrayObject<O, T>> getMedoids()
        {

            std::vector<dataset::BasicArrayObject<O, T>> medoids;

            for (size_t i = 0; i < clusters.size(); i++)
                medoids.push_back(clusters[i].getMedoid());

            return medoids;

        }

        void run(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                 std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nClusters,
                 size_t numberOfIterations, size_t seed)
        {

            std::vector<size_t> aux1, aux2;
            size_t ite = 0;

            initialize(dataset, nClusters, seed);

            while (ite < numberOfIterations)
            {

                aux1 = getMedoidsOID();

                assignment(dataset, df);

#ifdef ENABLE_DEBUG
                std::cout << "\n================= Iteration: " << ite << "=================" << std::endl;
                std::cout << *this;
                std::cout << "=========================================================================\n" << std::endl;
#endif

                reCenter(df);

                aux2 = getMedoidsOID();

                if (aux1 == aux2)
                    break;

                aux1.clear();
                aux2.clear();

                ite++;

            }

        }

        friend std::ostream& operator<<(std::ostream& os, const Kmedoids& kmedoids)
        {

            os << "Clusters: " << std::endl << std::endl;

            for (auto &cluster : kmedoids.clusters) {
                os << cluster << std::endl << std::endl;
            }

            return os;

        }

    };
}

#endif //GERVLIB_KMEDOIDS_H
