//
// Created by joaovictor on 15/08/23.
//

#include "Configure.h"
#include "Utils.h"
#include "Dataset.h"

using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;

int main(int argc, char **argv)
{

    auto data = std::make_unique<Dataset<size_t, double>>("../data/Dataset1.csv", " ");
    data->getElement(0).setOID(5);
    std::cout << *data << std::endl;

    data->clear();



    return 0;

}