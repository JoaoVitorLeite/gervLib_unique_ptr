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

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSeed(16);

    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 500, 8000, true, false, false, "tmp_unit_test11");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");

    std::vector<vptree::Node<size_t, double>*> nodes = vp->getLeafNodes();
    int sum = 0;
    for (auto node : nodes) {
        auto leaf = dynamic_cast<vptree::LeafNode<size_t, double>*>(node);
        if (leaf->getIndex() != nullptr)
            sum += leaf->getIndex()->getDataset()->getCardinality();
        else
            sum += leaf->getDataset()->getCardinality();
    }

    std::cout << "SUM = " << sum << std::endl;

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(i), 5, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 5, true);

        for(size_t j = 0; j < res1.size(); j++)
        {
            if (res1[j].getDistance() != res2[j].getDistance()) {
                //assert(res1[j].getDistance() == res2[j].getDistance());
                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
                throw std::runtime_error("Error");
            }
        }
//        std::cout << "LEAF = " << vp->getLeafNodeAccess() << std::endl;
    }

//    std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(47), 5, true);
//    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(6), 100, true);

//    for(size_t j = 0; j < res1.size(); j++)
//    {
//        if (res1[j].getDistance() != res2[j].getDistance()) {
//                //assert(res1[j].getDistance() == res2[j].getDistance());
//            std::cout << "Error index " << 0 << ": " << res1[j].getDistance() << " != " << res2[j].getDistance() << std::endl;
//            throw std::runtime_error("Error");
//        }
//    }

    gervLib::utils::deleteDirectory("tmp_unit_test11");
    gervLib::utils::deleteDirectory("tmp_unit_test12");
    return 0;


//    std::unique_ptr<Index<size_t, double>> vp2 = std::make_unique<vptree::VPTree<size_t, double>>();
//
//    std::unique_ptr<u_char[]> d = vp->serialize();
//
//    vp2->deserialize(std::move(d));
//
//    std::cout << std::endl << vp->isEqual(vp2) << std::endl;




    return 0;

}
