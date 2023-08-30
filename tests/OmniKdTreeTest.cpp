//
// Created by joaoleite on 8/29/23.
//

#include "OmniKdTree.h"
#include "LAESA.h"
#include <cassert>

using namespace gervLib::index;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int test1()
{

    auto pivot = gervLib::dataset::BasicArrayObject<size_t, double>(1, {1.0, 2.0, 3.0});
    std::unique_ptr<kdtree::Node<size_t, double>> node = std::make_unique<kdtree::Node<size_t, double>>();
    node->setBoundsSize(2);
    node->setBound(0, 0.0, 1.0);
    node->setBound(1, 0.14, 16.7);

    std::unique_ptr<u_char[]> serialized = node->serialize();

    std::unique_ptr<kdtree::Node<size_t, double>> deserialized = std::make_unique<kdtree::Node<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getBoundsSize() == 2);
    assert(deserialized->getBound(0).first == 0.0);
    assert(deserialized->getBound(0).second == 1.0);
    assert(deserialized->getBound(1).first == 0.14);
    assert(deserialized->getBound(1).second == 16.7);
    assert(deserialized->isEqual(node));

    return 0;
}

int test2()
{
    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    auto laesa = std::make_unique<LAESA<size_t, double>>(std::move(data), std::move(dist), std::move(pvt), 2, "tmp_unit_test18");
    laesa->saveIndex();
    std::unique_ptr<Index<size_t, double>> laesa2 = std::make_unique<LAESA<size_t, double>>("tmp_unit_test18", "");
    std::unique_ptr<kdtree::LeafNode<size_t, double>> leaf = std::make_unique<kdtree::LeafNode<size_t, double>>();
    leaf->setBoundsSize(2);
    leaf->setBound(0, 0.0, 1.0);
    leaf->setBound(1, 0.14, 16.7);
    leaf->setIndex(std::move(laesa));

    std::unique_ptr<u_char[]> serialized = leaf->serialize();
    std::unique_ptr<kdtree::LeafNode<size_t, double>> deserialized = std::make_unique<kdtree::LeafNode<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getIndex()->isEqual(laesa2));
    assert(deserialized->getBound(0).first == 0.0);
    assert(deserialized->getBound(0).second == 1.0);
    assert(deserialized->getBound(1).first == 0.14);
    assert(deserialized->getBound(1).second == 16.7);

    std::unique_ptr<kdtree::Node<size_t, double>> aux = std::move(leaf);
    assert(deserialized->isEqual(aux));

    gervLib::utils::deleteDirectory("tmp_unit_test18");
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

    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4096, false, true, true, "tmp_unit_test19");

    std::unique_ptr<u_char[]> serialized = omni->serialize();
    std::unique_ptr<Index<size_t, double>> deserialized = std::make_unique<omni::OmniKdTree<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(omni->isEqual(deserialized));

    gervLib::utils::deleteDirectory("tmp_unit_test19");
    return 0;
}

int test4()
{
    gervLib::configure::configure();
    std::cout << std::boolalpha;

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data2 = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            test = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ",");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
    pvt->setSeed(16);

    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test20");
    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test21");

    for(size_t i = 0; i < test->getCardinality(); i++)
    {
        std::vector<gervLib::query::ResultEntry<size_t>> res1 = omni->kNNIncremental(test->getElement(i), 100, true);
        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);

        for(size_t j = 0; j < res1.size(); j++)
            assert(res1[j].getDistance() == res2[j].getDistance());
    }

    gervLib::utils::deleteDirectory("tmp_unit_test20");
    gervLib::utils::deleteDirectory("tmp_unit_test21");
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