//
// Created by joaovictor on 09/08/23.
//

#include <cassert>
#include "PageManager.h"
#include "Page.h"

using namespace gervLib::memory;

int test1()
{

    Page<size_t> page = Page<size_t>(0, 10, "test.bin");

    std::unique_ptr<u_char[]> data = page.serialize();
    Page<size_t> page2 = Page<size_t>();
    page2.deserialize(std::move(data));

    assert(page.getStart() == page2.getStart());
    assert(page.getEnd() == page2.getEnd());
    assert(page.getPath() == page2.getPath());

    return 0;

}

int test2()
{

//    PageManager<size_t> pageManager = PageManager<size_t>("page", "tmp_unit_test", 5);
//    double data[] = {1.0, 2.0, 3.0, 4.0};
//    std::unique_ptr<u_char[]> data2(new u_char[4 * sizeof(double)]);
//    std::memcpy(data2.get(), data, 4 * sizeof(double));
//    pageManager.save(17, std::move(data2), 4 * sizeof(double));
//
//    double data4[] = {5.0, 6.0, 7.0, 8.0};
//    std::unique_ptr<u_char[]> data5(new u_char[4 * sizeof(double)]);
//    std::memcpy(data5.get(), data4, 4 * sizeof(double));
//    pageManager.save(23, std::move(data5), 4 * sizeof(double));
//
//    std::unique_ptr<u_char[]> t = pageManager.serialize();
//    auto pageManager2 = std::unique_ptr<PageManager<size_t>>();
//    pageManager2->deserialize(std::move(t));
//
//    std::unique_ptr<u_char[]> data3 = pageManager2->load(17);
//    assert(((double *) data3.get())[0] == 1.0);
//    assert(((double *) data3.get())[1] == 2.0);
//    assert(((double *) data3.get())[2] == 3.0);
//    assert(((double *) data3.get())[3] == 4.0);
//
//    std::unique_ptr<u_char[]> data6 = pageManager2->load(23);
//    assert(((double *) data6.get())[0] == 5.0);
//    assert(((double *) data6.get())[1] == 6.0);
//    assert(((double *) data6.get())[2] == 7.0);
//    assert(((double *) data6.get())[3] == 8.0);
//
//    assert(pageManager.isEqual(*pageManager2));

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

    std::unique_ptr<u_char[]> data7 = pageManager.load(17);
    assert(((double *) data7.get())[0] == 1.0);
    assert(((double *) data7.get())[1] == 2.0);
    assert(((double *) data7.get())[2] == 3.0);
    assert(((double *) data7.get())[3] == 4.0);

    std::unique_ptr<u_char[]> data8 = pageManager.load(23);
    assert(((double *) data8.get())[0] == 5.0);
    assert(((double *) data8.get())[1] == 6.0);
    assert(((double *) data8.get())[2] == 7.0);
    assert(((double *) data8.get())[3] == 8.0);

    return 0;
}

int main(int argc, char *argv[])
{

    int res = 0;
    res += test1();
    res += test2();

    return res == 0 ? 0 : 1;

}