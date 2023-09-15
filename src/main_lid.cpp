//
// Created by joaoleite on 9/15/23.
//

#include "argparse.hpp"
#include "Dataset.h"
#include "Pivot.h"
#include "Pivots.h"
#include "EuclideanDistance.h"
#include "Indexes.h"
#include "DatasetUtils.h"
#include "DistanceFactory.h"

using namespace gervLib::dataset;
using namespace gervLib::pivots;
using namespace gervLib::distance;
using namespace gervLib::index;

int main(int argc, char **argv)
{

    argparse::ArgumentParser program("main_lid");

    program.add_argument("-INDEX_FOLDER")
            .help("Path to the index folder")
            .required();

    program.add_argument("-DATASET")
            .help("Path to the dataset file")
            .required();

    program.add_argument("-DATASET_SEPARATOR")
            .help("Separator used in the dataset file")
            .required();

    program.add_argument("-DISTANCE_FUNCTION")
            .help("Distance function to be used")
            .default_value("EUCLIDEAN");

    program.add_argument("-K_MAX")
            .help("Maximum number of neighbors to be used in the LID calculation")
            .default_value("100");

    program.add_argument("-FILE_NAME")
        .help("Name of the file to be used in the output")
        .required();

    program.add_argument("-OUTPUT_PATH")
        .help("Path to the output directory")
        .required();

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>(program.get<std::string>("-DATASET"), program.get<std::string>("-DATASET_SEPARATOR"));
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distance_function;

    if(program.get<std::string>("-DISTANCE_FUNCTION") == "EUCLIDEAN")
    {

        distance_function = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    }
    else
    {

        throw std::invalid_argument("Invalid distance function");

    }

    std::unique_ptr<Index<size_t, double>> index = std::make_unique<gervLib::index::vptree::VPTree<size_t, double>>(program.get<std::string>("-INDEX_FOLDER"));

    calculateLID(data1, distance_function, index, std::stoul(program.get<std::string>("-K_MAX")), program.get<std::string>("-FILE_NAME"), program.get<std::string>("-OUTPUT_PATH"));

    return 0;
}