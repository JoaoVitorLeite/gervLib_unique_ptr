//
// Created by joaovictor on 15/08/23.
//

#include "Configure.h"
#include "Utils.h"
#include "Dataset.h"
#include "EditDistance.h"
#include "EuclideanDistance.h"
#include "RandomPivots.h"

using namespace gervLib::configure;
using namespace gervLib::utils;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::pivots;

int main(int argc, char **argv)
{

    Timer timer;
    timer.start();
    timer.stop();
    std::cout << timer.getElapsedTime() << std::endl;
    std::cout << timer.getElapsedTimeUser() << std::endl;

    return 0;

}