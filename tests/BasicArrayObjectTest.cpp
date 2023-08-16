//
// Created by joaoleite on 14/07/23.
//

#include "BasicArrayObject.h"
#include <cassert>

using namespace gervLib::dataset;

int test1()
{

    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>();
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>();

    return 0;

}

int test2()
{

    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, 8);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, 2);

    assert(b1.getOID() == 54);
    assert(b1.size() == 8);
    assert(b2.getOID() == 87);
    assert(b2.size() == 2);

    return 0;

}

int test3()
{

    auto* b1 = new BasicArrayObject<size_t, double>(54, 8);
    auto *b2 = new BasicArrayObject<size_t, std::vector<char>>(87, 2);

    delete b1;
    delete b2;

    return 0;

}

int test4()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    assert(b1.getOID() == 54);
    assert(b1.size() == 6);
    assert(b1.getData() == v);
    assert(b2.getOID() == 87);
    assert(b2.size() == 2);
    assert(b2.getData() == v2);

    return 0;

}

int test5()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    auto b1 = std::make_unique<BasicArrayObject<size_t, double>>(54, v);
    auto b2 = std::make_unique<BasicArrayObject<size_t, std::vector<char>>>(87, v2);

    b1.reset();
    assert(b1 == nullptr);
    assert(b2 != nullptr);

    return 0;

}

int test6()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1.setOID(1);
    b2.setOID(2);

    assert(b1.getOID() == 1);
    assert(b2.getOID() == 2);

    return 0;

}

int test7()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    std::vector<double> v3 = std::vector<double>{1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v4 = {{'a', 'b', 'c'}};
    b1.setData(v3);
    b2.setData(v4);

    assert(b1.getData() == v3);
    assert(b2.getData() == v4);

    return 0;

}

int test8()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    assert(b1[0] == 1.0);
    assert(b1[1] == 2.0);
    assert(b1[2] == 3.0);
    assert(b1[3] == 4.0);
    assert(b1[4] == 5.0);
    assert(b1[5] == 6.0);
    assert(b2[0] == v2[0]);
    assert(b2[1] == v2[1]);

    return 0;

}

int test9()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1[0] = 2.0;
    b1[1] = 3.0;
    b1[2] = 4.0;
    b1[3] = 5.0;
    b1[4] = 6.0;
    b1[5] = 7.0;
    b2[0] = v2[1];
    b2[1] = v2[0];

    assert(b1[0] == 2.0);
    assert(b1[1] == 3.0);
    assert(b1[2] == 4.0);
    assert(b1[3] == 5.0);
    assert(b1[4] == 6.0);
    assert(b1[5] == 7.0);
    assert(b2[0] == v2[1]);
    assert(b2[1] == v2[0]);

    return 0;

}

int test10()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    assert(b1.size() == 6);
    assert(b2.size() == 2);

    return 0;

}

int test11()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1.clear();
    b2.clear();

    assert(b1.size() == 0);
    assert(b2.size() == 0);

    return 0;

}

int test12()
{

    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}, {'d', 'e', 'f'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1.set(0, 2.0);
    b1.set(1, 3.0);
    b1.set(2, 4.0);
    b1.set(3, 5.0);
    b1.set(4, 6.0);
    b1.set(5, 7.0);
    b2.set(0, v2[1]);
    b2.set(1, v2[0]);

    assert(b1[0] == 2.0);
    assert(b1[1] == 3.0);
    assert(b1[2] == 4.0);
    assert(b1[3] == 5.0);
    assert(b1[4] == 6.0);
    assert(b1[5] == 7.0);
    assert(b2[0] == v2[1]);
    assert(b2[1] == v2[0]);

    return 0;

}

int test13()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1.set(4.0);
    b2.set(v2[0]);

    assert(b1[3] == 4.0);
    assert(b2[1] == v2[0]);

    return 0;

}

int test14()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    std::unique_ptr<BasicArrayObject<size_t, double>> b3 = b1.clone();
    std::unique_ptr<BasicArrayObject<size_t, std::vector<char>>> b4 = b2.clone();

    assert(b3->getOID() == b1.getOID());
    assert(b3->size() == b1.size());
    assert(b3->getData() == b1.getData());
    assert(b4->getOID() == b2.getOID());
    assert(b4->size() == b2.size());
    assert(b4->getData() == b2.getData());

    return 0;

}

int test15()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    assert(b1.isEqual(b3));
    assert(b2.isEqual(b4));

    return 0;

}

int test16()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b3[0] = 2.0;
    b4[0] = {'d', 'e', 'f'};

    assert(!b1.isEqual(b3));
    assert(!b2.isEqual(b4));

    return 0;

}

int test17()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    assert(b1 == b3);
    assert(b2 == b4);

    return 0;

}

int test18()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b3[0] = 2.0;
    b4[0] = {'d', 'e', 'f'};

    assert(b1 != b3);
    assert(b2 != b4);

    return 0;

}

int test19()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);
    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>(88, v2);

    assert(!(b1 < b3));
    assert(b2 < b4);
    assert(b1 <= b3);
    assert(b2 <= b4);
    assert(!(b1 > b3));
    assert(!(b2 > b4));
    assert(b1 >= b3);
    assert(!(b2 >= b4));

    return 0;

}

int test20()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    auto s1 = b1.serialize();
    auto s2 = b2.serialize();

    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>();
    b3.deserialize(std::move(s1));
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>();
    b4.deserialize(std::move(s2));

    assert(b1 == b3);
    assert(b2 == b4);

    return 0;

}

int test21()
{

    std::vector<double> v = {1.0, 2.0, 3.0};
    std::vector<std::vector<char>> v2 = {{'a', 'b', 'c'}};
    BasicArrayObject<size_t, double> b1 = BasicArrayObject<size_t, double>(54, v);
    BasicArrayObject<size_t, std::vector<char>> b2 = BasicArrayObject<size_t, std::vector<char>>(87, v2);

    b1.clear();
    b2.clear();

    auto s1 = b1.serialize();
    auto s2 = b2.serialize();

    BasicArrayObject<size_t, double> b3 = BasicArrayObject<size_t, double>();
    b3.deserialize(std::move(s1));
    BasicArrayObject<size_t, std::vector<char>> b4 = BasicArrayObject<size_t, std::vector<char>>();
    b4.deserialize(std::move(s2));

    assert(b1 == b3);
    assert(b2 == b4);

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
    res += test11();
    res += test12();
    res += test13();
    res += test14();
    res += test15();
    res += test16();
    res += test17();
    res += test18();
    res += test19();
    res += test20();
    res += test21();

    return res == 0 ? 0 : 1;

}