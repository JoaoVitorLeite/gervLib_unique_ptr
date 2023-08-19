//
// Created by joaovictor on 03/08/23.
//

#include "Query.h"
#include "Dataset.h"
#include <cassert>

using namespace gervLib::dataset;
using namespace gervLib::query;

//Test serialize and deserialize knn entry with type as size_t
int test1()
{
    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    ResultEntry<size_t> result = ResultEntry<size_t>();
    result.setElement(49);
    result.setDistance(0.17);

    std::unique_ptr<u_char[]> arr = result.serialize();
    ResultEntry<size_t> result2 = ResultEntry<size_t>();
    result2.deserialize(std::move(arr));

    assert(result == result2);
    assert(result2.getElement() == 49);
    assert(result2.getDistance() == 0.17);

    return 0;
}

//Test serialize and deserialize knn entry with type as size_t
int test2()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");
    ResultEntry<BasicArrayObject<size_t, double>> result = ResultEntry<BasicArrayObject<size_t, double>>();
    result.setElement(data[0]);
    result.setDistance(0.996);

    std::unique_ptr<u_char[]> arr = result.serialize();
    ResultEntry<BasicArrayObject<size_t, double>> result2 = ResultEntry<BasicArrayObject<size_t, double>>();
    result2.deserialize(std::move(arr));

    assert(result.isEqual(result2));
    assert(result2.getElement() == data[0]);
    assert(result2.getDistance() == 0.996);

    return 0;

}

//Test custom priority with result entry of type size_t
int test3()
{

    std::function<bool(const ResultEntry<size_t> &, const ResultEntry<size_t> &)> compare = [](const ResultEntry<size_t> &a, const ResultEntry<size_t> &b) {
        return a.getDistance() > b.getDistance();
    };

    CustomPriorityQueue<ResultEntry<size_t>> queue = CustomPriorityQueue<ResultEntry<size_t>>(compare);

    assert(queue.empty() == true);
    assert(queue.size() == 0);

    queue.push(ResultEntry<size_t>(1, 0.1));
    queue.push(ResultEntry<size_t>(2, 0.2));
    queue.push(ResultEntry<size_t>(3, 0.3));
    queue.push(ResultEntry<size_t>(4, 0.4));

    assert(queue.size() == 4);
    assert(queue.empty() == false);
    assert(queue.top().getElement() == 1);
    assert(queue.top().getDistance() == 0.1);

    queue.pop();

    assert(queue.size() == 3);

    queue.clear();

    assert(queue.size() == 0);

    return 0;

}

//Test custom priority with result entry of type BasicArrayObject<size_t, double>
int test4()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    std::function<bool(const ResultEntry<BasicArrayObject<size_t, double>> &, const ResultEntry<BasicArrayObject<size_t, double>> &)> compare = [](const ResultEntry<BasicArrayObject<size_t, double>> &a, const ResultEntry<BasicArrayObject<size_t, double>> &b) {
        return a.getDistance() > b.getDistance();
    };

    CustomPriorityQueue<ResultEntry<BasicArrayObject<size_t, double>>> queue = CustomPriorityQueue<ResultEntry<BasicArrayObject<size_t, double>>>(compare);

    assert(queue.empty() == true);
    assert(queue.size() == 0);

    queue.push(ResultEntry<BasicArrayObject<size_t, double>>(data[0], 0.1));
    queue.push(ResultEntry<BasicArrayObject<size_t, double>>(data[1], 0.2));
    queue.push(ResultEntry<BasicArrayObject<size_t, double>>(data[2], 0.3));
    queue.push(ResultEntry<BasicArrayObject<size_t, double>>(data[3], 0.4));

    assert(queue.size() == 4);
    assert(queue.empty() == false);
    assert(queue.top().getElement() == data[0]);
    assert(queue.top().getDistance() == 0.1);

    queue.pop();

    assert(queue.size() == 3);

    queue.clear();

    assert(queue.size() == 0);

    return 0;

}

//Test result of type size_t
int test5()
{

    Result<size_t> result = Result<size_t>();

    assert(result.size() == 0);
    assert(result.empty() == true);

    result.push(ResultEntry<size_t>(1, 0.1));
    result.push(ResultEntry<size_t>(2, 0.2));
    result.push(ResultEntry<size_t>(3, 0.3));
    result.push(ResultEntry<size_t>(4, 0.4));

    assert(result.size() == 4);
    assert(result.empty() == false);
    assert(result.top().getElement() == 4);
    assert(result.top().getDistance() == 0.4);

    result.pop();

    assert(result.size() == 3);

    result.clear();

    assert(result.size() == 0);

    return 0;

}

//Test result of type BasicArrayObject<size_t, double>
int test6()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    Result<BasicArrayObject<size_t, double>> result = Result<BasicArrayObject<size_t, double>>(true);

    assert(result.size() == 0);
    assert(result.empty() == true);

    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[0], 0.1));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[1], 0.2));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[2], 0.3));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[3], 0.4));

    assert(result.size() == 4);
    assert(result.empty() == false);
    assert(result.top().getElement() == data[0]);
    assert(result.top().getDistance() == 0.1);

    result.pop();

    assert(result.size() == 3);

    result.clear();

    assert(result.size() == 0);

    return 0;

}

//Test serialize and deserialize for result of type size_t
int test7()
{

    Result<size_t> result = Result<size_t>();

    result.push(ResultEntry<size_t>(1, 0.1));
    result.push(ResultEntry<size_t>(2, 0.2));
    result.push(ResultEntry<size_t>(3, 0.3));
    result.push(ResultEntry<size_t>(4, 0.4));

    std::unique_ptr<u_char[]> arr = result.serialize();
    Result<size_t> result2 = Result<size_t>();
    result2.deserialize(std::move(arr));

    assert(result.isEqual(result2));
    assert(result2.top().getElement() == 4);
    assert(result2.top().getDistance() == 0.4);

    return 0;

}

//Test serialize and deserialize for result of type BasicArrayObject<size_t, double>
int test8()
{

    Dataset<size_t, double> data = Dataset<size_t, double>("../../data/Dataset1.csv", " ");

    Result<BasicArrayObject<size_t, double>> result = Result<BasicArrayObject<size_t, double>>(true);

    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[0], 0.1));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[1], 0.2));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[2], 0.3));
    result.push(ResultEntry<BasicArrayObject<size_t, double>>(data[3], 0.4));

    std::unique_ptr<u_char[]> arr = result.serialize();
    Result<BasicArrayObject<size_t, double>> result2 = Result<BasicArrayObject<size_t, double>>();
    result2.deserialize(std::move(arr));

    assert(result.isEqual(result2));
    assert(result2.top().getElement() == data[0]);
    assert(result2.top().getDistance() == 0.1);

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

    return res == 0 ? 0 : 1;

}