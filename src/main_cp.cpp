//
// Created by joaovictor on 15/08/23.
//

#include "Configure.h"
#include "Utils.h"
#include "Dataset.h"
//#include "EditDistance.h"
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
#include "PMTree.h"
#include "SPBTree.h"
#include "LC.h"
#include "ISPivots.h"
#include "FFTPivots.h"
#include "DatasetUtils.h"

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

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ",");

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSampleSize(1.0);
    pvt->setSeed(449462);

    std::unique_ptr<Index<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 12, 50, 0, false, false, true, true);
    std::cout << *vp << std::endl;
//    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4096, 2, 2, 4, 2, false, true, true, true);
//    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4096, false, true, true);
//    std::unique_ptr<pmtree::PMTree<size_t, double>> pm = std::make_unique<pmtree::PMTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true);
//    std::unique_ptr<spbtree::SPBTree<size_t, double>> spb = std::make_unique<spbtree::SPBTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 50, 4096, false, true, true);
//    std::unique_ptr<lc::LC<size_t, double>> lc = std::make_unique<lc::LC<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 0, true, true, true, true);

//    vp->saveIndex();

//    splitByLID(data2, dist2, vp, 100, "cities", "../data");
//    splitTrainTest<size_t, double>("../data/cities_norm.csv", ",", 0.8, 176);

//    std::vector<gervLib::query::ResultEntry<size_t>> res1 = pm->kNNIncremental(test->getElement(0), 100, true, true);
//    std::vector<gervLib::query::ResultEntry<size_t>> res2 = vp->kNNIncremental(test->getElement(0), 100, true, true);
//    std::cout << *lc << std::endl;
//
//    std::unique_ptr<u_char[]> serialized = lc->serialize();
//    std::unique_ptr<Index<size_t, double>> lc2 = std::make_unique<lc::LC<size_t, double>>();
//    lc2->deserialize(std::move(serialized));
//
//    std::cout << lc->isEqual(lc2) << std::endl;

//    std::unique_ptr<Index<size_t, double>> spb = std::make_unique<spbtree::SPBTree<size_t, double, mpz_class>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 5, 4096, false, false, true);

//    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
//    std::unique_ptr<spbtree::SPBTree<size_t, double>> spb2 = std::make_unique<spbtree::SPBTree<size_t, double>>();
//    std::unique_ptr<u_char[]> serialized = spb->serialize();
//    spb2->deserialize(std::move(serialized));

//    std::cout << *spb << "\n\n";

//    for(size_t i = 0; i < test->getCardinality(); i++)
//    {
//
//        std::vector<gervLib::query::ResultEntry<size_t>> res1 = lc2->kNNIncremental(test->getElement(i), 100, true, true);
//        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true, true);
//
//        for(size_t j = 0; j < res1.size(); j++)
//        {
//            if (res1[j].getDistance() != res2[j].getDistance()) {
//                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
//                //throw std::runtime_error("Error");
//            }
//
////            std::cout << spb->getPrunning() << std::endl;
//
//        }
//    }

//    std::vector<gervLib::query::ResultEntry<size_t>> res1 = lc->kNNIncremental(test->getElement(0), 5, true, true);
//    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(0), 5, true);
//
//    for(size_t j = 0; j < res1.size(); j++)
//    {
//        std::cout << res2[j] << "\t" << res1[j] << std::endl;
//    }

//    size_t id = 2;
//    std::vector<gervLib::query::ResultEntry<size_t>> res1 = spb->kNNIncremental(test->getElement(id), 10, true);
//    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(id), 10, true);
//
//    for(size_t j = 0; j < res1.size(); j++)
//    {
//        std::cout << res1[j] << "\t" << res2[j] << std::endl;
//        if (res1[j].getDistance() != res2[j].getDistance()) {
//                //assert(res1[j].getDistance() == res2[j].getDistance());
//            std::cout << "Error: " << res1[j].getDistance() << " != " << res2[j].getDistance() << std::endl;
//            throw std::runtime_error("Error");
//        }
//    }
//
//    std::vector<gervLib::query::ResultEntry<size_t>> res3 = spb->kNNIncremental(test->getElement(id), 10, true);
//    std::vector<gervLib::query::ResultEntry<size_t>> res4 = sc->kNN(test->getElement(id), 10, true);
//
//    for(size_t j = 0; j < res3.size(); j++)
//    {
//        std::cout << res3[j] << "\t" << res4[j] << std::endl;
//        if (res3[j].getDistance() != res4[j].getDistance()) {
//            //assert(res1[j].getDistance() == res2[j].getDistance());
//            std::cout << "Error: " << res3[j].getDistance() << " != " << res4[j].getDistance() << std::endl;
//            throw std::runtime_error("Error");
//        }
//    }


//    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test11");
//    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
//
//    for(size_t i = 0; i < test->getCardinality(); i++)
//    {
//        std::vector<gervLib::query::ResultEntry<size_t>> res1 = omni->kNNIncremental(test->getElement(i), 100, true);
//        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
//
//        for(size_t j = 0; j < res1.size(); j++)
//        {
//            if (res1[j].getDistance() != res2[j].getDistance()) {
//                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
//                //throw std::runtime_error("Error");
//            }
////            std::cout << omni->getPrunning() << std::endl;
//        }
//    }

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
