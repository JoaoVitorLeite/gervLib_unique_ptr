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
#include <cassert>

using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;
using namespace gervLib::memory;

int main(int argc, char **argv)
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    PageManager<size_t> pageManager = PageManager<size_t>("page", "tmp_unit_test", 5);
    double data[] = {1.0, 2.0, 3.0, 4.0};
    std::unique_ptr<u_char[]> data2(new u_char[4 * sizeof(double)]);
    std::memcpy(data2.get(), data, 4 * sizeof(double));
    pageManager.save(17, std::move(data2), 4 * sizeof(double));

    double data4[] = {5.0, 6.0, 7.0, 8.0};
    std::unique_ptr<u_char[]> data5(new u_char[4 * sizeof(double)]);
    std::memcpy(data5.get(), data4, 4 * sizeof(double));
    pageManager.save(23, std::move(data5), 4 * sizeof(double));

    std::unique_ptr<u_char[]> data3 = pageManager.load(17);
    assert(((double *) data3.get())[0] == 1.0);
    assert(((double *) data3.get())[1] == 2.0);
    assert(((double *) data3.get())[2] == 3.0);
    assert(((double *) data3.get())[3] == 4.0);

    std::unique_ptr<u_char[]> data6 = pageManager.load(23);
    assert(((double *) data6.get())[0] == 5.0);
    assert(((double *) data6.get())[1] == 6.0);
    assert(((double *) data6.get())[2] == 7.0);
    assert(((double *) data6.get())[3] == 8.0);

    std::unique_ptr<u_char[]> s = pageManager.serialize();
    PageManager<size_t> pageManager2 = PageManager<size_t>();
    pageManager2.deserialize(std::move(s));

    assert(pageManager.isEqual(pageManager2));


    return 0;

}