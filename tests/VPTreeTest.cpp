//
// Created by joaovictor on 06/08/23.
//

#include "VPTree.h"
#include "LAESA.h"
#include <cassert>

using namespace gervLib::index;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int test1()
{

    auto pivot = gervLib::dataset::BasicArrayObject<size_t, double>(1, {1.0, 2.0, 3.0});
    std::unique_ptr<vptree::Node<size_t, double>> node = std::make_unique<vptree::Node<size_t, double>>();
    node->setPivot(std::make_unique<BasicArrayObject<size_t, double>>(pivot));
    node->setMu(1.0);
    node->setCoverage(4.0);

    std::unique_ptr<u_char[]> serialized = node->serialize();

    std::unique_ptr<vptree::Node<size_t, double>> deserialized = std::make_unique<vptree::Node<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    std::cout << *deserialized << std::endl;

    assert(deserialized->getPivot()->getOID() == 1);
    assert(deserialized->getPivot()->getData() == std::vector<double>({1.0, 2.0, 3.0}));
    assert(deserialized->getMu() == 1.0);
    assert(deserialized->getCoverage() == 4.0);
    assert(deserialized->getPivot()->isEqual(pivot));
    assert(deserialized->isEqual(node));

    return 0;
}

int test2()
{
    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    auto laesa = std::make_unique<LAESA<size_t, double>>(std::move(data), std::move(dist), std::move(pvt), 2, "tmp_unit_test8");
    laesa->saveIndex();
    std::unique_ptr<Index<size_t, double>> laesa2 = std::make_unique<LAESA<size_t, double>>("tmp_unit_test8", "");
    std::unique_ptr<vptree::LeafNode<size_t, double>> leaf = std::make_unique<vptree::LeafNode<size_t, double>>();
    leaf->setCoverage(7.39);
    leaf->setIndex(std::move(laesa));

    std::unique_ptr<u_char[]> serialized = leaf->serialize();
    std::unique_ptr<vptree::LeafNode<size_t, double>> deserialized = std::make_unique<vptree::LeafNode<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getMu() == 0.0);
    assert(deserialized->getCoverage() == 7.39);
    assert(deserialized->getIndex()->isEqual(laesa2));

    std::unique_ptr<vptree::Node<size_t, double>> aux = std::move(leaf);
    assert(deserialized->isEqual(aux));

    gervLib::utils::deleteDirectory("tmp_unit_test8");
    return 0;
}

int test3()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " "),
                                             data2 = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " "),
                                             test = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
                                                                        dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();

    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 10, 0, false, true, true, true, "tmp_unit_test9");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test10");


    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(i), 5, true, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 5, true, true);

        for(size_t j = 0; j < res1.size(); j++)
            assert(res1[j].getDistance() == res2[j].getDistance());

    }

    gervLib::utils::deleteDirectory("tmp_unit_test9");
    gervLib::utils::deleteDirectory("tmp_unit_test10");
    return 0;

}

int test4()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();

    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4900, false, true, true, true, "tmp_unit_test11");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(i), 100, true, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true, true);

        for(size_t j = 0; j < res1.size(); j++)
            assert(res1[j].getDistance() == res2[j].getDistance());

    }

    gervLib::utils::deleteDirectory("tmp_unit_test11");
    gervLib::utils::deleteDirectory("tmp_unit_test12");
    return 0;

}

int test5()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();

    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4900, false, true, true, true, "tmp_unit_test13");

    std::unique_ptr<u_char[]> serialized = vp->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<vptree::VPTree<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(vp->isEqual(deserialized));

    gervLib::utils::deleteDirectory("tmp_unit_test13");
    return 0;

}

int test6()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();

    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4900, false, true, true, true, "tmp_unit_test31");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test32");

    std::unique_ptr<u_char[]> serialized = vp->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<vptree::VPTree<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(vp->isEqual(deserialized));

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = deserialized->kNNIncremental(test->getElement(i), 100, true, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true, true);

        for(size_t j = 0; j < res1.size(); j++)
            assert(res1[j].getDistance() == res2[j].getDistance());

    }

    gervLib::utils::deleteDirectory("tmp_unit_test31");
    gervLib::utils::deleteDirectory("tmp_unit_test32");
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
    res += test4();
    res += test5();

    return res == 0 ? 0 : 1;

}