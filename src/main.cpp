//
// Created by joaovictor on 15/08/23.
//

#include "Configure.h"
#include "Utils.h"
#include "Dataset.h"
#include "EditDistance.h"

using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;
using namespace gervLib::distance;

int main(int argc, char **argv)
{

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, std::vector<char>>>> df = std::make_unique<EditDistance<BasicArrayObject<size_t, std::vector<char>>>>();
    Dataset<size_t, std::vector<char>> data = Dataset<size_t, std::vector<char>>("../data/names.csv", " ");

    double d = df->operator()(data[0], data[1]);
    std::cout << "\nDIST = " << d << std::endl;
    std::cout << data[0] << std::endl;
    std::cout << data[1] << std::endl;



    return 0;

}