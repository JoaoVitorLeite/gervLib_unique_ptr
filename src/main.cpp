//
// Created by joaovictor on 15/08/23.
//

#include "Configure.h"
#include "Utils.h"
#include "Dataset.h"
#include "EditDistance.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"
#include "KmedoidsPivots.h"
#include "PageManager.h"
#include "SequentialScan.h"
#include <cassert>

using namespace gervLib::index;
using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;
using namespace gervLib::memory;
using namespace gervLib::query;

int main(int argc, char **argv)
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    BasicArrayObject<unsigned long, double> obj = data->getElement(0);
    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data), std::move(distanceFunction), "tmp_unit_test3");

    std::vector<ResultEntry<size_t>> results = sequentialScan->kNN(obj, 5, true);

    std::unique_ptr<u_char[]> serialized = sequentialScan->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<SequentialScan<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    std::vector<ResultEntry<size_t>> results2 = deserialized->kNN(obj, 5, true);

    assert(results[0].getElement() == results2[0].getElement());
    assert(results[1].getElement() == results2[1].getElement());
    assert(results[2].getElement() == results2[2].getElement());
    assert(results[3].getElement() == results2[3].getElement());
    assert(results[4].getElement() == results2[4].getElement());

    assert(results[0].getDistance() == results2[0].getDistance());
    assert(results[1].getDistance() == results2[1].getDistance());
    assert(results[2].getDistance() == results2[2].getDistance());
    assert(results[3].getDistance() == results2[3].getDistance());
    assert(results[4].getDistance() == results2[4].getDistance());

    assert(sequentialScan->isEqual(deserialized));

    gervLib::utils::deleteDirectory("tmp_unit_test3");


    return 0;

}