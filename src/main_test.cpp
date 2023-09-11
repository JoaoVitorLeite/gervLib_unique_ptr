//
// Created by joaoleite on 9/11/23.
//

#include <iostream>
#include "Indexes.h"
#include "Pivots.h"
#include "Query.h"

/* KeyWords:
 * -INDEX:
 * -DATASET_TRAIN:
 *      -DATASET_TRAIN_CARDINALITY:
 *      -DATASET_TRAIN_DIMENSIONALITY:
 *      -DATASET_TRAIN_SEPARATOR:
 * -DATASET_TEST:
 *      -DATASET_TEST_CARDINALITY:
 *      -DATASET_TEST_DIMENSIONALITY:
 *      -DATASET_TEST_SEPARATOR:
 * -DISTANCE_FUNCTION: Default Euclidean
 * -PIVOT_TYPE: Default RandomPivots
 *      -PIVOT_SAMPLE_SIZE:
 *      -NUM_PIVOTS:
 * -SEED: Default 0
 * -K_MAX: Default 100
 * -NUM_QUERY: Default 100
 * -NUM_PER_LEAF: Default 1 % of the dataset train cardinality
 * -NUM_BINS: Default 200
 * -PAGE_SIZE: Default 0
 * -NUM_QUERIES_PER_FILE: Default 10000
 * -STORE_PIVOT_LEAF: Default true
 * -STORE_DIRECTORY_NODE: Default false
 * -STORE_LEAF_NODE: Default true
 * -USE_LAESA: Default true
 * */

using namespace gervLib::dataset;
using namespace gervLib::index;
using namespace gervLib::pivots;
using namespace gervLib::distance;
using namespace gervLib::query;

int main(int argc, char** argv)
{

    if((argc-1) % 2 != 0)
    {

        throw std::invalid_argument("Invalid number of arguments !_!");

    }
    else
    {

        std::map<std::string, std::string> args;
        size_t seed, k_max, num_query, num_per_leaf, num_bins, page_size, num_queries_per_file, num_pivots;
        bool store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa;
        std::variant<size_t, double> pivot_sample_size;

        std::unique_ptr<Dataset<size_t, double>> dataset_train, dataset_test;
        std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distance_function;
        std::unique_ptr<Pivot<size_t, double>> pivot;
        std::unique_ptr<Index<size_t, double>> index;

        for(int i = 1; i < argc; i += 2) {

            std::string key = argv[i];

            std::string value = argv[i + 1];

            for (char & x : key)
                x = std::toupper(x);

            args[key] = value;

        }

        if (args.find("-INDEX") == args.end() || args.find("-DATASET_TRAIN") == args.end() || args.find("-DATASET_TEST") == args.end() || args.find("-PIVOT_TYPE") == args.end())
        {

            throw std::invalid_argument("Index, dataset train, dataset test or pivot type not specified");

        }

        //Dataset train
        if (args.find("-DATASET_TRAIN_CARDINALITY") != args.end() && args.find("-DATASET_TRAIN_DIMENSIONALITY") != args.end())
        {
            dataset_train = std::make_unique<Dataset<size_t, double>>(args["-DATASET_TRAIN"], std::stoul(args["-DATASET_TRAIN_CARDINALITY"]), std::stoul(args["-DATASET_TRAIN_DIMENSIONALITY"]));
        }
        else if (args.find("-DATASET_TRAIN_SEPARATOR") != args.end())
        {
            dataset_train = std::make_unique<Dataset<size_t, double>>(args["-DATASET_TRAIN"], args["-DATASET_TRAIN_SEPARATOR"]);
        }
        else
            throw std::invalid_argument("Dataset train cardinality, dataset train dimensionality or dataset train separator not specified");

        //Dataset test
        if (args.find("-DATASET_TEST_CARDINALITY") != args.end() && args.find("-DATASET_TEST_DIMENSIONALITY") != args.end())
        {
            dataset_test = std::make_unique<Dataset<size_t, double>>(args["-DATASET_TEST"], std::stoul(args["-DATASET_TEST_CARDINALITY"]), std::stoul(args["-DATASET_TEST_DIMENSIONALITY"]));
        }
        else if (args.find("-DATASET_TEST_SEPARATOR") != args.end())
        {
            dataset_test = std::make_unique<Dataset<size_t, double>>(args["-DATASET_TEST"], args["-DATASET_TEST_SEPARATOR"]);
        }
        else
            throw std::invalid_argument("Dataset test cardinality, dataset test dimensionality or dataset test separator not specified");

        //Distance Function
        if (args.find("-DISTANCE_FUNCTION") == args.end())
        {

            distance_function = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

        }
        else
        {

            if(args["-DISTANCE_FUNCTION"] == "EUCLIDEAN")
            {

                distance_function = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

            }
            else
            {

                throw std::invalid_argument("Invalid distance function");

            }

        }

        //Seed
        if (args.find("-SEED") == args.end())
        {

            seed = 0;

        }
        else
        {

            seed = std::stoul(args["-SEED"]);

        }

        //K_max
        if (args.find("-K_MAX") == args.end())
        {

            k_max = 100;

        }
        else
        {

            k_max = std::stoul(args["-K_MAX"]);

        }

        //Num_query
        if (args.find("-NUM_QUERY") == args.end())
        {

            num_query = 100;

        }
        else
        {

            num_query = std::stoul(args["-NUM_QUERY"]);

        }

        //Num_per_leaf
        if (args.find("-NUM_PER_LEAF") == args.end())
        {
            num_per_leaf = std::ceil(0.01 * (double) dataset_train->getCardinality());
        }
        else
        {

            num_per_leaf = std::stoul(args["-NUM_PER_LEAF"]);

        }

        //Num_bins
        if (args.find("-NUM_BINS") == args.end())
        {

            num_bins = 200;

        }
        else
        {

            num_bins = std::stoul(args["-NUM_BINS"]);

        }

        //Page_size
        if (args.find("-PAGE_SIZE") == args.end())
        {

            page_size = 0;

        }
        else
        {

            page_size = std::stoul(args["-PAGE_SIZE"]);

        }

        //Num_queries_per_file
        if (args.find("-NUM_QUERIES_PER_FILE") == args.end())
        {

            num_queries_per_file = 10000;

        }
        else
        {

            num_queries_per_file = std::stoul(args["-NUM_QUERIES_PER_FILE"]);

        }

        //Store_pivot_leaf
        if (args.find("-STORE_PIVOT_LEAF") == args.end())
        {

            store_pivot_leaf = true;

        }
        else
        {

            if(args["-STORE_PIVOT_LEAF"] == "TRUE")
            {

                store_pivot_leaf = true;

            }
            else if(args["-STORE_PIVOT_LEAF"] == "FALSE")
            {

                store_pivot_leaf = false;

            }
            else
            {

                throw std::invalid_argument("Invalid store pivot leaf");

            }

        }

        //Store_directory_node
        if (args.find("-STORE_DIRECTORY_NODE") == args.end())
        {

            store_directory_node = false;

        }
        else
        {

            if(args["-STORE_DIRECTORY_NODE"] == "TRUE")
            {

                store_directory_node = true;

            }
            else if(args["-STORE_DIRECTORY_NODE"] == "FALSE")
            {

                store_directory_node = false;

            }
            else
            {

                throw std::invalid_argument("Invalid store directory node");

            }

        }

        //Store_leaf_node
        if (args.find("-STORE_LEAF_NODE") == args.end())
        {

            store_leaf_node = true;

        }
        else
        {

            if(args["-STORE_LEAF_NODE"] == "TRUE")
            {

                store_leaf_node = true;

            }
            else if(args["-STORE_LEAF_NODE"] == "FALSE")
            {

                store_leaf_node = false;

            }
            else
            {

                throw std::invalid_argument("Invalid store leaf node");

            }

        }

        //Use_laesa
        if (args.find("-USE_LAESA") == args.end())
        {

            use_laesa = true;

        }
        else
        {

            if(args["-USE_LAESA"] == "TRUE")
            {

                use_laesa = true;

            }
            else if(args["-USE_LAESA"] == "FALSE")
            {

                use_laesa = false;

            }
            else
            {

                throw std::invalid_argument("Invalid use laesa");

            }

        }

        //Num_pivots
        if (args.find("-NUM_PIVOTS") == args.end())
        {

            throw std::invalid_argument("Num pivots not specified");

        }
        else
        {

            num_pivots = std::stoul(args["-NUM_PIVOTS"]);

        }

        //Pivot_sample_size
        if (args.find("-PIVOT_SAMPLE_SIZE") == args.end())
        {

            throw std::invalid_argument("Pivot sample size not specified");

        }
        else
        {

            pivot_sample_size = std::stod(args["-PIVOT_SAMPLE_SIZE"]);

        }

        //Pivot_type
        if (args.find("-PIVOT_TYPE") == args.end())
        {
            pivot = std::make_unique<RandomPivots<size_t, double>>();
        }
        else
        {

            if (args["-PIVOT_TYPE"] == "RANDOM")
            {
                pivot = std::make_unique<RandomPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "CONVEX")
            {
                pivot = std::make_unique<ConvexPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "FFT")
            {
                pivot = std::make_unique<FFTPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "BPP")
            {
                pivot = std::make_unique<BPPPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "HFI")
            {
                pivot = std::make_unique<HFIPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "IS")
            {
                pivot = std::make_unique<ISPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "KMEDOIDS")
            {
                pivot = std::make_unique<KmedoidsPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "MAXSEPARATED")
            {
                pivot = std::make_unique<MaxSeparatedPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "MAXVARIANCE")
            {
                pivot = std::make_unique<MaxVariancePivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "PCA")
            {
                pivot = std::make_unique<PCAPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "SELECTION")
            {
                pivot = std::make_unique<SelectionPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "WDR")
            {
                pivot = std::make_unique<WDRPivots<size_t, double>>();
            }
            else if (args["-PIVOT_TYPE"] == "SSS")
            {
                pivot = std::make_unique<SSSPivots<size_t, double>>();
            }
            else
            {
                throw std::invalid_argument("Invalid pivot type");
            }

        }

        pivot->setSeed(seed);
        pivot->setSampleSize(pivot_sample_size);

        //Index
        if (args["-INDEX"] == "VPTREE")
        {
            index = std::make_unique<gervLib::index::vptree::VPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa);
        }
        else if (args["-INDEX"] == "MVPTREE")
        {
            index = std::make_unique<gervLib::index::mvptree::MVPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, 2, 2, 4, 2, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa);
        }
        else if (args["-INDEX"] == "OMNIKDTREE")
        {
            index = std::make_unique<gervLib::index::omni::OmniKdTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa);
        }
        else if (args["-INDEX"] == "PMTREE")
        {
            index = std::make_unique<gervLib::index::pmtree::PMTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa);
        }
        else if (args["-INDEX"] == "SPBTREE")
        {
            index = std::make_unique<gervLib::index::spbtree::SPBTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, num_bins, page_size, store_directory_node, store_leaf_node, use_laesa);
        }

        //Queries
        size_t max_queries = std::min(dataset_test->getCardinality(), num_query);
        std::vector<ResultEntry<size_t>> results;

        for (size_t i = 0; i < dataset_test->getCardinality(); i++)
        {
            index->kNNIncremental(dataset_test->getElement(i), k_max, false, true, results);

            if ((i + 1) % num_queries_per_file == 0)
                index->generateIndexFiles(false, true);
        }

    }

    return 0;

}