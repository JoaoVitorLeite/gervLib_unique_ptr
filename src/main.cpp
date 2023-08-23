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
#include "VPTree.h"

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

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<RandomPivots<size_t, double>>();
    std::unique_ptr<vptree::VPTree<size_t, double>> vp =
            std::make_unique<vptree::VPTree<size_t, double>>(std::move(data), std::move(distanceFunction),
                    std::move(pivots), 2, 12, 0, false, true, true);

    std::unique_ptr<Index<size_t, double>> vp2 = std::make_unique<vptree::VPTree<size_t, double>>();

    std::unique_ptr<u_char[]> d = vp->serialize();

    vp2->deserialize(std::move(d));

    std::cout << std::endl << vp->isEqual(vp2) << std::endl;

//    vp2->deserialize(std::move(d));

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



    return 0;

}
