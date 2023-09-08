//
// Created by joaovictor on 08/09/23.
//

#include "Configure.h"
#include "Dataset.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"
#include "LAESA.h"
#include "LC.h"

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

    std::unique_ptr<lc::LC<size_t, double>> lc = std::make_unique<lc::LC<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 0, true, true, true, true, "tmp_unit_test28");

    std::unique_ptr<u_char[]> serialized = lc->serialize();
    std::unique_ptr<Index<size_t, double>> lc2 = std::make_unique<lc::LC<size_t, double>>();
    lc2->deserialize(std::move(serialized));

    assert(lc->isEqual(lc2));

    gervLib::utils::deleteDirectory("tmp_unit_test28");
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

    std::unique_ptr<lc::LC<size_t, double>> lc = std::make_unique<lc::LC<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 0, true, true, true, true, "tmp_unit_test29");

    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test30");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = lc->kNNIncremental(test->getElement(i), 100, true, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true, true);

        for(size_t j = 0; j < res1.size(); j++)
        {
            assert(res1[j].getDistance() == res2[j].getDistance());
        }
    }

    gervLib::utils::deleteDirectory("tmp_unit_test29");
    gervLib::utils::deleteDirectory("tmp_unit_test30");
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