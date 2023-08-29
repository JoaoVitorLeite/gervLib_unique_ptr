//
// Created by joaovictor on 26/08/23.
//

#include "Configure.h"
#include "MVPTree.h"
#include "Dataset.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"
#include "LAESA.h"

using namespace gervLib::pivots;
using namespace gervLib::distance;
using namespace gervLib::index;
using namespace gervLib::dataset;

int test1()
{
    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<mvptree::Node<size_t, double>> node = std::make_unique<mvptree::Node<size_t, double>>(2, 4, 2);
    node->setPivot(0, data->getElement(0));
    node->setPivot(1, data->getElement(5));
    node->setNodeID(8);
    node->setMemoryStatus(MEMORY_STATUS::IN_MEMORY);
    node->setSplit(0, 0, 0.27);
    node->setSplit(1, 0, 0.187);
    node->setSplit(1, 1, 0.468);

    assert(node->getNumberOfPivots() == 2);
    assert(node->getNumberOfSplits() == 2);
    assert(node->getNodeID() == 8);
    assert(node->getPivot(0).isEqual(data->getElement(0)));
    assert(node->getPivot(1).isEqual(data->getElement(5)));

    std::unique_ptr<u_char[]> serialized = node->serialize();
    std::unique_ptr<mvptree::Node<size_t, double>> deserialized = std::make_unique<mvptree::Node<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getNumberOfPivots() == 2);
    assert(deserialized->getNumberOfSplits() == 2);
    assert(deserialized->getNodeID() == 8);
    assert(deserialized->getPivot(0).isEqual(data->getElement(0)));
    assert(deserialized->getPivot(1).isEqual(data->getElement(5)));
    assert(deserialized->getSplit(0, 0) == 0.27);
    assert(deserialized->getSplit(1, 0) == 0.187);
    assert(deserialized->getSplit(1, 1) == 0.468);
    assert(deserialized->getMemoryStatus() == MEMORY_STATUS::IN_MEMORY);

    return 0;
}

int test2()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    auto obj1 = data->getElement(0), obj2 = data->getElement(5);
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    auto laesa = std::make_unique<LAESA < size_t, double>>(std::move(data), std::move(dist), std::move(pvt), 2, "tmp_unit_test13");
    laesa->saveIndex();
    std::unique_ptr<Index<size_t, double>> laesa2 = std::make_unique<LAESA<size_t, double>>("tmp_unit_test13", "");
    std::unique_ptr<mvptree::LeafNode<size_t, double>> leaf = std::make_unique<mvptree::LeafNode < size_t, double>>();
    leaf->setNumberOfPivots(2);
    leaf->setPivot(0, obj1);
    leaf->setPivot(1, obj2);
    leaf->setNodeID(8);
    leaf->setMemoryStatus(MEMORY_STATUS::IN_MEMORY);
    leaf->setNumberOfSplits(2, 2);
    leaf->setSplit(0, 0, 0.27);
    leaf->setSplit(1, 0, 0.187);
    leaf->setSplit(1, 1, 0.468);
    leaf->setIndex(std::move(laesa));

    std::unique_ptr<u_char[]> serialized = leaf->serialize();
    std::unique_ptr<mvptree::LeafNode<size_t, double>> deserialized = std::make_unique<mvptree::LeafNode<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getNumberOfPivots() == 2);
    assert(deserialized->getNumberOfSplits() == 2);
    assert(deserialized->getNodeID() == 8);
    assert(deserialized->getPivot(0).isEqual(obj1));
    assert(deserialized->getPivot(1).isEqual(obj2));
    assert(deserialized->getSplit(0, 0) == 0.27);
    assert(deserialized->getSplit(1, 0) == 0.187);
    assert(deserialized->getSplit(1, 1) == 0.468);
    assert(deserialized->getMemoryStatus() == MEMORY_STATUS::IN_MEMORY);
    assert(deserialized->getIndex()->isEqual(laesa2));

    gervLib::utils::deleteDirectory("tmp_unit_test13");
    return 0;

}

int test3()
{
    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();

    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2,
                                                                                                               500, 4096, 2, 2, 4, 2, false, true, true, false, "tmp_unit_test14");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test15");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = mvp->kNNIncremental(test->getElement(i), 100, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);

        for(size_t j = 0; j < res1.size(); j++)
            assert(res1[j].getDistance() == res2[j].getDistance());

    }

    gervLib::utils::deleteDirectory("tmp_unit_test14");
    gervLib::utils::deleteDirectory("tmp_unit_test15");
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

    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2,
                                                                                                               500, 4096, 2, 2, 4, 2, false, true, true, false, "tmp_unit_test16");

    std::unique_ptr<u_char []> serialized = mvp->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<mvptree::MVPTree<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(mvp->isEqual(deserialized));

    gervLib::utils::deleteDirectory("tmp_unit_test16");
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

    return res == 0 ? 0 : 1;

}