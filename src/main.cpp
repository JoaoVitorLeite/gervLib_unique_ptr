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
#include "MVPTree.h"
#include "OmniKdTree.h"
//#include "KdTree.h"

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

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSeed(16);

    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test11");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = omni->kNNIncremental(test->getElement(i), 100, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);

        for(size_t j = 0; j < res1.size(); j++)
        {
            if (res1[j].getDistance() != res2[j].getDistance()) {
                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
                //throw std::runtime_error("Error");
            }
//            std::cout << omni->getPrunning() << std::endl;
        }
    }

//    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 4096);
//    std::unique_ptr<u_char[]> serialized = mvp->serialize();
//    std::unique_ptr<Index<size_t, double>> mvp2 = std::make_unique<mvptree::MVPTree<size_t, double>>();
//    mvp2->deserialize(std::move(serialized));
//
//    std::cout << mvp->isEqual(mvp2) << std::endl;

//    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
//            data2 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
//            test = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ",");
//    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
//            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
//
//    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
//
//    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2,
//            500, 4096, 2, 2, 4, 2, false, true, true, false, "tmp_unit_test11");
//    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
//
//    for(size_t i = 0; i < test->getCardinality(); i++)
//    {
//        std::vector<gervLib::query::ResultEntry<size_t>> res1 = mvp->kNNIncremental(test->getElement(i), 100, true);
//        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
//
//        for(size_t j = 0; j < res1.size(); j++)
//            assert(res1[j].getDistance() == res2[j].getDistance());
//
//        std::cout << mvp->getPrunning() << std::endl;
//
//    }
//
//    gervLib::utils::deleteDirectory("tmp_unit_test11");
//    gervLib::utils::deleteDirectory("tmp_unit_test12");


//    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test11");
//    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
//
//    for(size_t i = 0; i < test->getCardinality(); i++)
//    {
//        std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(i), 100, true);
//        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
//
//        for(size_t j = 0; j < res1.size(); j++)
//        {
//            if (res1[j].getDistance() != res2[j].getDistance()) {
//                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
//                throw std::runtime_error("Error");
//            }
//            std::cout << vp->getPrunning() << std::endl;
//        }
//    }

////    std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(47), 5, true);
////    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(6), 100, true);
//
////    for(size_t j = 0; j < res1.size(); j++)
////    {
////        if (res1[j].getDistance() != res2[j].getDistance()) {
////                //assert(res1[j].getDistance() == res2[j].getDistance());
////            std::cout << "Error index " << 0 << ": " << res1[j].getDistance() << " != " << res2[j].getDistance() << std::endl;
////            throw std::runtime_error("Error");
////        }
////    }
//
//    gervLib::utils::deleteDirectory("tmp_unit_test11");
//    gervLib::utils::deleteDirectory("tmp_unit_test12");
//
    return 0;


}
