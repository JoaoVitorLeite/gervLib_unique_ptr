//
// Created by joaovictor on 01/08/23.
//

#include <cassert>
#include "Dataset.h"
#include "EuclideanDistance.h"
#include "SelectionPivots.h"

using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int test1()
{

    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<SelectionPivots<size_t, double>>();
    pvt->setSeed(837);
    pvt->operator()(data, dist, 4);
    std::unique_ptr<u_char[]> arr = pvt->serialize();
    std::unique_ptr<Pivot<unsigned long, double>> pvt2 = std::make_unique<SelectionPivots<size_t, double>>();
    pvt2->deserialize(std::move(arr));
    assert(pvt->isEqual(pvt2));

    return 0;

}

int test2()
{

#ifndef ENABLE_SRAND
    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<SelectionPivots<size_t, double>>();
    pvt->setSeed(837);
    pvt->operator()(data, dist, 4);
    assert(pvt->getPivot(0).getOID() == 6);
    assert(pvt->getPivot(1).getOID() == 8);
    assert(pvt->getPivot(2).getOID() == 11);
    assert(pvt->getPivot(3).getOID() == 18);
#else
    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<SelectionPivots<size_t, double>>();
    pvt->setSeed(837);
    pvt->operator()(data, dist, 4);
    assert(pvt->getPivot(0).getOID() == 7);
    assert(pvt->getPivot(1).getOID() == 10);
    assert(pvt->getPivot(2).getOID() == 13);
    assert(pvt->getPivot(3).getOID() == 16);
#endif

    return 0;

}

int test3()
{

#ifndef ENABLE_SRAND
    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<SelectionPivots<size_t, double>>();
    pvt->setSeed(837);
    pvt->setNumberOfDropPivots(2);
    pvt->operator()(data, dist, 2);
    assert(pvt->getPivot(0).getOID() == 11);
    assert(pvt->getPivot(1).getOID() == 18);
#else
    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<SelectionPivots<size_t, double>>();
    pvt->setSeed(837);
    pvt->setNumberOfDropPivots(2);
    pvt->operator()(data, dist, 2);
    assert(pvt->getPivot(0).getOID() == 13);
    assert(pvt->getPivot(1).getOID() == 16);
#endif

    return 0;

}

int main(int argc, char *argv[])
{

    int res = 0;
    res += test1();
    res += test2();
    res += test3();

    return res == 0 ? 0 : 1;

}