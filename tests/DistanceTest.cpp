//
// Created by joaoleite on 16/07/23.
//

#include "DistanceFunction.h"
#include "EuclideanDistance.h"
//#include "EditDistance.h"
#include "Dataset.h"
#include "DistanceFactory.h"
#include <cassert>
#include <cmath>

using namespace gervLib::dataset;
using namespace gervLib::distance;

int test1()
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    return 0;

}

int test2()
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    assert(df->getDistanceType() == EUCLIDEAN);
    assert(df->operator()(data[0], data[0]) == 0.0);
    assert(df->operator()(data[0], data[1]) == sqrt(37));
    assert(df->operator()(data[0], data[3]) == sqrt(5));
    assert(df->operator()(data[5], data[12]) == df->operator()(data[12], data[5]));

    return 0;

}

int test3()
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    std::unique_ptr<u_char[]> serialized = df->serialize();
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    df2->deserialize(std::move(serialized));

    assert(df2->getDistanceCount() == 0);
    assert(df->isEqual(df2));
    assert(df->operator()(data[5], data[12]) == df2->operator()(data[12], data[5]));

    return 0;

}

int test4()
{
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df2 = DistanceFactory<BasicArrayObject<size_t, double>>::createDistanceFunction(df->getDistanceType());

    assert(df->isEqual(df2));

    return 0;

}

int main(int argc, char *argv[])
{

    int res = 0;
    res += test1();
    res += test2();
    res += test3();
    res += test4();

    return res == 0 ? 0 : 1;

}