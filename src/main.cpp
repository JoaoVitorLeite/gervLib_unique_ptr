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

using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int main(int argc, char **argv)
{

    auto data = std::make_unique<Dataset<size_t, double>>("../data/Dataset1.csv", 20, 2);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<KmedoidsPivots<size_t, double>>();
    pvt->operator()(data, dist, 4);

    return 0;

}