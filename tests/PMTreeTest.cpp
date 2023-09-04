//
// Created by joaoleite on 9/1/23.
//

#include "Configure.h"
#include "Dataset.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"
#include "LAESA.h"
#include "PMTree.h"

using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;
using namespace gervLib::index;

int test1()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSeed(16);

    std::unique_ptr<pmtree::PMTree<size_t, double>> pm = std::make_unique<pmtree::PMTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt),
                                                                                                          2, 50, 8000, false, true, true, "tmp_unit_test22");

    std::unique_ptr<u_char[]> serialized = pm->serialize();
    std::unique_ptr<Index<size_t, double>> pm2 = std::make_unique<pmtree::PMTree<size_t, double>>();
    pm2->deserialize(std::move(serialized));

    assert(pm->isEqual(pm2));

    gervLib::utils::deleteDirectory("tmp_unit_test22");
    return 0;
}

int test2()
{

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSeed(16);

    std::unique_ptr<pmtree::PMTree<size_t, double>> pm = std::make_unique<pmtree::PMTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt),
                                                                                                          2, 50, 8000, false, true, true, "tmp_unit_test23");

    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test24");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = pm->kNNIncremental(test->getElement(i), 100, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);

        for(size_t j = 0; j < res1.size(); j++)
        {
            assert(res1[j].getDistance() == res2[j].getDistance());
        }
    }

    gervLib::utils::deleteDirectory("tmp_unit_test23");
    gervLib::utils::deleteDirectory("tmp_unit_test24");
    return 0;

}

int main(int argc, char *argv[])
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    int res = 0;
    res += test1();
    res += test2();

    return res == 0 ? 0 : 1;

}