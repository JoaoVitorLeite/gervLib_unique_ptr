//
// Created by joaoleite on 16/07/23.
//

#include "DistanceFunction.h"
#include "EuclideanDistance.h"
#include "EditDistance.h"
#include "Dataset.h"
#include "DistanceFactory.h"
#include <cassert>
#include <cmath>

using namespace gervLib::dataset;
using namespace gervLib::distance;

int test1()
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto df2 = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();

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

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();
    Dataset<size_t, std::vector<char>> data = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    assert(df->getDistanceType() == LEVENSHTEIN);
    assert(df->operator()(data[0], data[0]) == 0.0);
    assert(df->operator()(data[0], data[1]) == 8.0);
    assert(df->operator()(data[3], data[10]) == 5.0);
    assert(df->operator()(data[5], data[12]) == df->operator()(data[12], data[5]));

    return 0;

}

int test4()
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

int test5()
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();
    Dataset<size_t, std::vector<char>> data = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    double d;
    d = df->operator()(data[0], data[0]);
    assert(d == 0.0);
    d = df->operator()(data[0], data[1]);
    assert(d == 8.0);
    d = df->operator()(data[3], data[10]);
    assert(d == 5.0);
    d = df->operator()(data[5], data[12]);
    assert(d == df->operator()(data[12], data[5]));

    std::unique_ptr<u_char[]> serialized = df->serialize();
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df2 = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();
    df2->deserialize(std::move(serialized));

    assert(df2->getDistanceCount() == 5);
    assert(df->isEqual(df2));
    assert(df->operator()(data[5], data[12]) == df2->operator()(data[12], data[5]));

    return 0;

}

int test6()
{
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df2 = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df3 = DistanceFactory<BasicArrayObject<size_t, double>>::createDistanceFunction(df->getDistanceType());
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df4 = DistanceFactory<BasicArrayObject<size_t, std::vector<char>>>::createDistanceFunction(df2->getDistanceType());

    assert(df->isEqual(df3));
    assert(df2->isEqual(df4));

    return 0;

}

int main(int argc, char *argv[])
{

    int res = 0;
    res += test1();
    res += test2();
    res += test3();
    res += test4();
    res += test5();
    res += test6();

    return res == 0 ? 0 : 1;

}