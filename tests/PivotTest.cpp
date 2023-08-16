//
// Created by joaovictor on 01/08/23.
//

#include <cassert>
#include "Dataset.h"
#include "DistanceFunction.h"
#include "EuclideanDistance.h"
#include "Pivot.h"

using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int test1()
{

    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<Pivot<size_t, double>>();

    assert(pvt->getNumberOfPivots() == 0);
    assert(pvt->getPivotType() == PIVOT_TYPE::RANDOM);
    pvt->setSeed(167);
    assert(pvt->getSeed() == 167);

    pvt->getPivots()->setData(std::vector<BasicArrayObject<size_t, double>>(2));
    pvt->setPivot(0, data->getElement(0));
    pvt->setPivot(1, data->getElement(5));

    data->clear();

    assert(pvt->getNumberOfPivots() == 2);

    return 0;

}

int test2()
{

    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<Pivot<size_t, double>>();

    assert(pvt->getNumberOfPivots() == 0);
    assert(pvt->getPivotType() == PIVOT_TYPE::RANDOM);
    pvt->setSeed(167);
    assert(pvt->getSeed() == 167);

    pvt->getPivots()->setData(std::vector<BasicArrayObject<size_t, double>>(2));
    pvt->setPivot(0, data->getElement(0));
    pvt->setPivot(1, data->getElement(5));

    data->clear();

    assert(pvt->getNumberOfPivots() == 2);

    std::unique_ptr<u_char[]> buffer = pvt->serialize();
    auto pvt2 = std::make_unique<Pivot<size_t, double>>();
    pvt2->deserialize(std::move(buffer));

    assert(pvt2->getNumberOfPivots() == 2);
    assert(pvt2->getPivotType() == PIVOT_TYPE::RANDOM);
    assert(pvt2->getSeed() == 167);
    assert(pvt2->getPivot(0).getOID() == 0);
    assert(pvt2->getPivot(1).getOID() == 5);
    assert(pvt2->isEqual(pvt));

    return 0;

}

int test3()
{

    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<Pivot<size_t, double>>();

    assert(pvt->getNumberOfPivots() == 0);
    assert(pvt->getPivotType() == PIVOT_TYPE::RANDOM);
    pvt->setSeed(167);
    assert(pvt->getSeed() == 167);

    pvt->getPivots()->setData(std::vector<BasicArrayObject<size_t, double>>(2));
    pvt->setPivot(0, data->getElement(0));
    pvt->setPivot(1, data->getElement(5));

    data->clear();

    assert(pvt->getNumberOfPivots() == 2);

    pvt->clear();

    std::unique_ptr<u_char[]> buffer = pvt->serialize();
    auto pvt2 = std::make_unique<Pivot<size_t, double>>();
    pvt2->deserialize(std::move(buffer));

    assert(pvt2->getNumberOfPivots() == 0);
    assert(pvt2->getPivotType() == PIVOT_TYPE::RANDOM);
    assert(pvt2->getSeed() == 167);
    assert(pvt2->isEqual(pvt));

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