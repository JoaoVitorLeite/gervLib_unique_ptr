//
// Created by joaoleite on 9/1/23.
//

#include "Configure.h"
#include "Dataset.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"
#include "LAESA.h"
#include "PMTree.h"

using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;
using namespace gervLib::index;

int test1()
{

    auto pivot = gervLib::dataset::BasicArrayObject<size_t, double>(1, {1.0, 2.0, 3.0});
    std::unique_ptr<pmtree::Node<size_t, double>> node = std::make_unique<pmtree::Node<size_t, double>>();
    node->setPivot(std::make_unique<BasicArrayObject<size_t, double>>(pivot));
    node->setCoverage(4.0);

    std::unique_ptr<u_char[]> serialized = node->serialize();

    std::unique_ptr<pmtree::Node<size_t, double>> deserialized = std::make_unique<pmtree::Node<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getPivot()->getOID() == 1);
    assert(deserialized->getPivot()->getData() == std::vector<double>({1.0, 2.0, 3.0}));
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
    std::unique_ptr<pmtree::LeafNode<size_t, double>> leaf = std::make_unique<pmtree::LeafNode<size_t, double>>();
    leaf->setCoverage(7.39);
    leaf->setIndex(std::move(laesa));
    leaf->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

    std::unique_ptr<u_char[]> serialized = leaf->serialize();
    std::unique_ptr<pmtree::LeafNode<size_t, double>> deserialized = std::make_unique<pmtree::LeafNode<size_t, double>>();
    deserialized->deserialize(std::move(serialized));

    assert(deserialized->getCoverage() == 7.39);
    assert(deserialized->getIndex()->isEqual(laesa2));

    std::unique_ptr<pmtree::Node<size_t, double>> aux = std::move(leaf);
    assert(deserialized->isEqual(aux));

    gervLib::utils::deleteDirectory("tmp_unit_test8");
    return 0;
}

int main(int argc, char *argv[])
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    int res = 0;
    res += test1();
    res += test2();

    return res == 0 ? 0 : 1;

}