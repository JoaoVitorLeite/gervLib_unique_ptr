//
// Created by joaoleite on 14/07/23.
//

#include "Dataset.h"
#include "BasicArrayObject.h"
#include "DatasetWrapper.h"
#include <cassert>

using namespace gervLib::dataset;

int test1()
{

    Dataset<size_t, double> data = Dataset<size_t, double>();
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>();
    return 0;

}

int test2()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    return 0;

}

int test3()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, double> data2 = Dataset<size_t, double>(data);

    assert(data == data2);
    assert(data.getSeed() == data2.getSeed());
    assert(data.getCardinality() == data2.getCardinality());
    assert(data.getDimensionality() == data2.getDimensionality());
    assert(data.getPath() == data2.getPath());

    return 0;

}

int test4()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, double> data2 = data;

    assert(data == data2);
    assert(data.getSeed() == data2.getSeed());
    assert(data.getCardinality() == data2.getCardinality());
    assert(data.getDimensionality() == data2.getDimensionality());
    assert(data.getPath() == data2.getPath());

    return 0;

}

int test5()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    assert(data.getCardinality() == 20);
    assert(data2.getCardinality() == 2472);

    return 0;

}

int test6()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    assert(data.getDimensionality() == 2);
    assert(data2.getDimensionality() == 2);

    return 0;

}

int test7()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    assert(data.getSeed() == 0);
    assert(data2.getSeed() == 0);

    return 0;

}

int test8()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    assert(data.getPath() == "../../data/Dataset1.csv");
    assert(data2.getPath() == "../../data/names.csv");

    return 0;

}

int test9()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, std::vector<char>> data2 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");

    data.setDimensionality(10);
    data2.setDimensionality(100);

    assert(data.getDimensionality() == 10);
    assert(data2.getDimensionality() == 100);

    return 0;

}

int test10()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    data.setSeed(10);

    assert(data.getSeed() == 10);

    return 0;

}

int test11()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    data.setPath("../../data/Dataset2.csv");

    assert(data.getPath() == "../../data/Dataset2.csv");

    return 0;

}

int test12()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, double> data2 = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    Dataset<size_t, double> data3 = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    data3.setDimensionality(80);

    assert(data == data2);
    assert(data != data3);

    return 0;

}

int test13()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    std::vector<double> v = {13.0, 5.0};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(0, v);
    BasicArrayObject<size_t, double> b2 = BasicArrayObject<size_t, double>(1, v);

    assert(data.contains(b1));
    assert(!data.contains(b2));

    return 0;

}

int test14()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(0, {13.0, 5.0});
    BasicArrayObject<size_t, double> b2 = BasicArrayObject<size_t, double>(1, {19.0, 4.0});
    data.erase(b1);

    assert(!data.contains(b1));
    assert(data.contains(b2));

    return 0;

}

int test15()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(0, {13.0, 5.0});
    BasicArrayObject<size_t, double> b2 = BasicArrayObject<size_t, double>(1, {19.0, 4.0});
    data.erase(0);

    assert(!data.contains(b1));
    assert(data.contains(b2));

    return 0;

}

int test16()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    std::unique_ptr<u_char[]> v = data.serialize();
    Dataset<size_t, double> data2 = Dataset<size_t, double>();
    data2.deserialize(std::move(v));

    Dataset<size_t, std::vector<char>> data3 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");
    std::unique_ptr<u_char[]> v2 = data3.serialize();
    Dataset<size_t, std::vector<char>> data4 = Dataset<size_t, std::vector<char>>();
    data4.deserialize(std::move(v2));

    assert(data == data2);
    assert(data3 == data4);

    return 0;

}

int test17()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    data.clear();
    std::unique_ptr<u_char[]> v = data.serialize();
    Dataset<size_t, double> data2 = Dataset<size_t, double>();
    data2.deserialize(std::move(v));

    Dataset<size_t, std::vector<char>> data3 = Dataset<size_t, std::vector<char>>("../../data/names.csv", " ");
    data.clear();
    std::unique_ptr<u_char[]> v2 = data3.serialize();
    Dataset<size_t, std::vector<char>> data4 = Dataset<size_t, std::vector<char>>();
    data4.deserialize(std::move(v2));

    assert(data == data2);
    assert(data3 == data4);

    return 0;

}

int test18()
{

    auto data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    data->erase(0);
    data->erase(data->operator[](1));

    assert(data->getCardinality() == 18);
    assert(data->operator[](0).getOID() == 1);
    assert(data->operator[](1).getOID() == 3);

    DatasetWrapper<size_t, double> data2 = DatasetWrapper<size_t, double>(data);

    assert(data->operator[](0).getOID() == 0);
    assert(data->operator[](1).getOID() == 1);

    data2.resetIndex(data);

    assert(data->operator[](0).getOID() == 1);
    assert(data->operator[](1).getOID() == 3);

    return 0;

}

int main(int argc, char *argv[])
{

    int res = 0;
    res += test1();
    res += test2();
    res += test3();
    res += test4();
    res += test5();
    res += test6();
    res += test7();
    res += test8();
    res += test9();
    res += test10();
    res += test11();
    res += test12();
    res += test13();
    res += test14();
    res += test15();
    res += test16();
    res += test17();
    res += test18();

    return res == 0 ? 0 : 1;

}

