//
// Created by joaovictor on 04/08/23.
//

#include "SequentialScan.h"
#include <cassert>

using namespace gervLib::index;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::query;

//Test sequential scan query for element 0
int test1()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    BasicArrayObject<unsigned long, double> obj = data->getElement(0);
    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data), std::move(distanceFunction), "tmp_unit_test1");

    std::vector<ResultEntry<size_t>> results = sequentialScan->kNN(obj, 5, true, true);

    assert(results[0].getElement() == 0);
    assert(results[1].getElement() == 3);
    assert(results[2].getElement() == 4);
    assert(results[3].getElement() == 16);
    assert(results[4].getElement() == 5);

    assert(results[0].getDistance() == 0.0);
    assert(results[1].getDistance() == sqrt(5.0));
    assert(results[2].getDistance() == sqrt(17.0));
    assert(results[3].getDistance() == sqrt(26.0));
    assert(results[4].getDistance() == sqrt(29.0));

    gervLib::utils::deleteDirectory("tmp_unit_test1");

    return 0;

}

//Test sequential scan query for element 6
int test2()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    BasicArrayObject<unsigned long, double> obj = data->getElement(6);
    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data), std::move(distanceFunction), "tmp_unit_test2");

    std::vector<ResultEntry<size_t>> results = sequentialScan->kNN(obj, 5, true, true);

    assert(results[0].getElement() == 6);
    assert(results[1].getElement() == 13);
    assert(results[2].getElement() == 9);
    assert(results[3].getElement() == 5);
    assert(results[4].getElement() == 7);
    assert(results[5].getElement() == 11);
    assert(results[6].getElement() == 3);

    assert(results[0].getDistance() == 0.0);
    assert(results[1].getDistance() == 2.0);
    assert(results[2].getDistance() == sqrt(80.0));
    assert(results[3].getDistance() == sqrt(82.0));
    assert(results[4].getDistance() == sqrt(148.0));
    assert(results[5].getDistance() == sqrt(149.0));
    assert(results[6].getDistance() == sqrt(160.0));

    gervLib::utils::deleteDirectory("tmp_unit_test2");

    return 0;

}

//Test sequential scan serialization
int test3()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    BasicArrayObject<unsigned long, double> obj = data->getElement(0);
    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data), std::move(distanceFunction), "tmp_unit_test3");

    std::vector<ResultEntry<size_t>> results = sequentialScan->kNN(obj, 5, true, true);

    std::unique_ptr<u_char[]> serialized = sequentialScan->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<SequentialScan<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    std::vector<ResultEntry<size_t>> results2 = deserialized->kNN(obj, 5, true, true);

    assert(results[0].getElement() == results2[0].getElement());
    assert(results[1].getElement() == results2[1].getElement());
    assert(results[2].getElement() == results2[2].getElement());
    assert(results[3].getElement() == results2[3].getElement());
    assert(results[4].getElement() == results2[4].getElement());

    assert(results[0].getDistance() == 0.0);
    assert(results[1].getDistance() == sqrt(5.0));
    assert(results[2].getDistance() == sqrt(17.0));
    assert(results[3].getDistance() == sqrt(26.0));
    assert(results[4].getDistance() == sqrt(29.0));

    assert(results[0].getDistance() == results2[0].getDistance());
    assert(results[1].getDistance() == results2[1].getDistance());
    assert(results[2].getDistance() == results2[2].getDistance());
    assert(results[3].getDistance() == results2[3].getDistance());
    assert(results[4].getDistance() == results2[4].getDistance());

    assert(results2[0].getDistance() == 0.0);
    assert(results2[1].getDistance() == sqrt(5.0));
    assert(results2[2].getDistance() == sqrt(17.0));
    assert(results2[3].getDistance() == sqrt(26.0));
    assert(results2[4].getDistance() == sqrt(29.0));

    assert(sequentialScan->isEqual(deserialized));

    gervLib::utils::deleteDirectory("tmp_unit_test3");

    return 0;

}

int main(int argc, char *argv[])
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    int res = 0;

    res += test1();
    res += test2();
    res += test3();

    return res == 0 ? 0 : 1;

}