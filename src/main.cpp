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
#include "LAESA.h"

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

//    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../data/Dataset1.csv", " ");
//    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
//    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<KmedoidsPivots<size_t, double>>();
//    std::unique_ptr<Index<size_t, double>> index = std::make_unique<LAESA<size_t, double>>(std::move(data), std::move(distanceFunction), std::move(pivots), 2, "tmp_laesa");
//
//    std::cout << *index << std::endl;
//
//    index->saveIndex();
//
////    index->clear();
//
//    std::unique_ptr<Index<size_t, double>> index2 = std::make_unique<LAESA<size_t, double>>("tmp_laesa", "");
//
//    std::cout << *index2 << std::endl;
//    std::cout << index->isEqual(index2) << std::endl;

    // Create a unique_ptr and allocate a new object
    std::unique_ptr<int> original_ptr = std::make_unique<int>(42);

    // Clone the object managed by original_ptr
    std::unique_ptr<int> cloned_ptr = std::make_unique<int>(*original_ptr);

    cloned_ptr.reset();
    std::cout << *original_ptr << std::endl;
    original_ptr.reset();

    return 0;

}