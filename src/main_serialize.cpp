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

    argparse::ArgumentParser program("main_serialize_vptree_maxvar");

    program.add_argument("-DATASET")
            .help("Path to the dataset file")
            .required();

    program.add_argument("-DATASET_SEPARATOR")
            .help("Separator used in the dataset file")
            .required();

    program.add_argument("-DISTANCE_FUNCTION")
            .help("Distance function to be used")
            .default_value("EUCLIDEAN");

    program.add_argument("-PIVOT_SAMPLE_SIZE")
            .help("Number of objects to be used in the pivot sample")
            .required();

    program.add_argument("-NUM_PIVOTS")
            .help("Number of pivots to be used")
            .required();

    program.add_argument("-SEED")
            .help("Seed for random number generation")
            .default_value("42");

    program.add_argument("-K_MAX")
            .help("Maximum number of neighbors to be used in the LID calculation")
            .default_value("100");

    program.add_argument("-NUM_PER_LEAF")
            .help("Number of objects per leaf");

    program.add_argument("-PAGE_SIZE")
            .help("Page size for the index")
            .default_value("0");

    program.add_argument("-STORE_PIVOT_LEAF")
            .help("Store the pivot in leaf")
            .default_value(true)
            .implicit_value(true);

    program.add_argument("-STORE_DIRECTORY_NODE")
            .help("Store the directory node")
            .default_value(false)
            .implicit_value(false);

    program.add_argument("-STORE_LEAF_NODE")
            .help("Store the leaf node")
            .default_value(false)
            .implicit_value(false);

    program.add_argument("-USE_LAESA")
            .help("Use LAESA in leaf")
            .default_value(false)
            .implicit_value(false);

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
    std::variant<size_t, double> pivot_sample_size;

    if(program.get<std::string>("-DISTANCE_FUNCTION") == "EUCLIDEAN")
    {

        distance_function = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    }
    else
    {

        throw std::invalid_argument("Invalid distance function");

    }

    std::filesystem::path path(gervLib::configure::baseOutputPath);
    path /= "VPTREE_MAX_VAR_SERIALIZE";

    std::unique_ptr<Pivot<size_t, double>> pivot = std::make_unique<MaxVariancePivots<size_t, double>>();
    pivot->setSampleSize(std::stod(program.get<std::string>("-PIVOT_SAMPLE_SIZE")));
    pivot->setSeed(std::stoul(program.get<std::string>("-SEED")));

    std::unique_ptr<Index<size_t, double>> index = 
            std::make_unique<gervLib::index::vptree::VPTree<size_t, double>>(std::move(data1), std::move(distance_function), std::move(pivot),
                    std::stoul(program.get<std::string>("-NUM_PIVOTS")), std::stoul(program.get<std::string>("-NUM_PER_LEAF")),
                    std::stoul(program.get<std::string>("-NUM_PER_LEAF")), program.get<bool>("-STORE_PIVOT_LEAF"),
                    program.get<bool>("-STORE_DIRECTORY_NODE"), program.get<bool>("-STORE_LEAF_NODE"),
                    program.get<bool>("-USE_LAESA"), path.string());

    index->saveIndex();
    
    return 0;
}