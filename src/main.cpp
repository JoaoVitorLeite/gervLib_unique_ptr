////
//// Created by joaovictor on 15/08/23.
////
//
//#include "Configure.h"
//#include "Dataset.h"
//#include "EuclideanDistance.h"
//#include "Indexes.h"
//#include "Pivots.h"
//#include "argparse.hpp"
//
//using namespace gervLib::index;
//using namespace gervLib::configure;
//using namespace gervLib::utils;
//using namespace gervLib::dataset;
//using namespace gervLib::distance;
//using namespace gervLib::pivots;
//using namespace gervLib::memory;
//using namespace gervLib::query;
//
//int main(int argc, char **argv)
//{
//
//    gervLib::configure::configure();
//    std::cout << std::boolalpha;
//
//    argparse::ArgumentParser program("gervLib");
//
//    program.add_argument("-op", "--operation")
//            .help("Operation to be performed (SERIALIZATION, CALCULATE_LID, QUERY)")
//            .required();
//
//    program.add_argument("--index")
//            .help("Index to be used")
//            .required();
//
//    program.add_argument("--dataset_train")
//            .help("Dataset to be used for training");
//
//    program.add_argument("--dataset_train_cardinality")
//            .help("Cardinality of the dataset to be used for training");
//
//    program.add_argument("--dataset_train_dimensionality")
//            .help("Dimensionality of the dataset to be used for training");
//
//    program.add_argument("--dataset_train_separator")
//            .help("Separator of the dataset to be used for training")
//            .default_value("");
//
//    program.add_argument("--dataset_test")
//            .help("Dataset to be used for testing");
//
//    program.add_argument("--dataset_test_cardinality")
//            .help("Cardinality of the dataset to be used for testing");
//
//    program.add_argument("--dataset_test_dimensionality")
//            .help("Dimensionality of the dataset to be used for testing");
//
//    program.add_argument("--dataset_test_separator")
//            .help("Separator of the dataset to be used for testing")
//            .default_value("");
//
//    program.add_argument("--distance_function")
//            .help("Distance function to be used")
//            .default_value("EUCLIDEAN");
//
//    program.add_argument("--pivot_type")
//            .help("Pivot type to be used")
//            .required();
//
//    program.add_argument("--pivot_sample_size")
//            .help("Pivot sample size to be used")
//            .default_value("1.0");
//
//    program.add_argument("--num_pivots")
//            .help("Number of pivots to be used");
//
//    //TODO - Add more arguments for pivot types (e.g. --pivot_candidates_size)
//
//    program.add_argument("--seed")
//            .help("Seed to be used")
//            .default_value("42");
//
//    program.add_argument("-k", "--knn")
//            .help("Number of nearest neighbors to be used")
//            .default_value("100");
//
//    program.add_argument("--num_queries")
//            .help("Number of queries to be executed");
//
//    program.add_argument("--num_per_leaf")
//            .help("Number of elements per leaf to be used")
//            .required();
//
//    program.add_argument("--num_bins")
//            .help("Number of bins to be used")
//            .default_value("256");
//
//    program.add_argument("--page_size")
//            .help("Page size to be used")
//            .default_value("0");
//
//    program.add_argument("--num_queries_per_file")
//            .help("Number of queries per file")
//            .default_value("10000");
//
//    program.add_argument("--store_pivot_leaf")
//            .help("Store pivot leaf")
//            .default_value("true");
//
//    program.add_argument("--store_directory_node")
//            .help("Store directory node")
//            .default_value("false");
//
//    program.add_argument("--store_leaf_node")
//            .help("Store leaf node")
//            .default_value("true");
//
//    program.add_argument("--use_laesa")
//            .help("Use LAESA in leaf nodes")
//            .default_value("true");
//
//    program.add_argument("--index_folder")
//            .help("Path to the index folder");
//
//    program.add_argument("--lid_file_name")
//            .help("Name of the file to be used in the output");
//
//    program.add_argument("--lid_output_path")
//            .help("Path to the output directory");
//
//    program.add_argument("--serialize_folder_name")
//            .help("Name of the folder to be used in the output")
//            .default_value("");
//
//    program.add_argument("--serialize_output_path")
//            .help("Path to the output directory")
//            .default_value("");
//
//    try {
//        program.parse_args(argc, argv);
//    }
//    catch (const std::runtime_error& err) {
//        std::cerr << err.what() << std::endl;
//        std::cerr << program;
//        std::exit(1);
//    }
//
//    std::unique_ptr<Dataset<size_t, double>> dataset_train, dataset_test;
//    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distance_function;
//    std::unique_ptr<Pivot<size_t, double>> pivot;
//    std::unique_ptr<Index<size_t, double>> index;
//    size_t k, seed, num_queries, num_per_leaf, num_bins, page_size, num_queries_per_file, num_pivots;
//
//    if (program.get<std::string>("-op") == "SERIALIZATION")
//    {
//
//        if (program.get("--dataset_train_separator") != "")
//        {
//            dataset_train = std::make_unique<Dataset<size_t, double>>(program.get<std::string>("--dataset_train"), program.get<std::string>("--dataset_train_separator"));
//        }
//        else
//        {
//            dataset_train = std::make_unique<Dataset<size_t, double>>(program.get<std::string>("--dataset_train"), std::stoul(program.get<std::string>("--dataset_train_cardinality")), std::stoul(program.get<std::string>("--dataset_train_dimensionality")));
//        }
//
//        if (program.get("--distance_function") == "EUCLIDEAN")
//        {
//            distance_function = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
//        }
//        else
//            throw std::invalid_argument("Invalid distance function");
//
//        PIVOT_TYPE pivot_type = gervLib::pivots::STR2PIVOT_TYPE[program.get<std::string>("--pivot_type")];
//        pivot = PivotFactory<size_t, double>::createPivot(pivot_type);
//
//        seed = std::stoul(program.get<std::string>("--seed"));
//        pivot->setSeed(seed);
//
//        if (program.get("--pivot_sample_size").find(".") != std::string::npos)
//        {
//            pivot->setSampleSize(std::stod(program.get<std::string>("--pivot_sample_size")));
//        }
//        else
//        {
//            pivot->setSampleSize(std::stoul(program.get<std::string>("--pivot_sample_size")));
//        }
//
//        //More pivot arguments
//
//        std::filesystem::path path;
//
//        if (program.get("--serialize_output_path") != "")
//            path = baseOutputPath;
//        else
//            path = program.get<std::string>("--serialize_output_path");
//
//        if (program.get("--serialize_folder_name") != "")
//            path /= program.get<std::string>("--serialize_folder_name");
//        else
//            path /= program.get("--index") + "_" + program.get("--pivot_type") + "_SERIALIZE";
//
//        gervLib::utils::deleteDirectory(path.string());
//
//        num_pivots = std::stoul(program.get<std::string>("--num_pivots"));
//        num_per_leaf = std::stoul(program.get<std::string>("--num_per_leaf"));
//        num_bins = std::stoul(program.get<std::string>("--num_bins"));
//        page_size = std::stoul(program.get<std::string>("--page_size"));
//        bool store_pivot_leaf = program.get<std::string>("--store_pivot_leaf") == "true";
//        bool store_directory_node = program.get<std::string>("--store_directory_node") == "true";
//        bool store_leaf_node = program.get<std::string>("--store_leaf_node") == "true";
//        bool use_laesa = program.get<std::string>("--use_laesa") == "true";
//
//        if (program.get("--index") == "VPTREE")
//        {
//            index = std::make_unique<gervLib::index::vptree::VPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function),std::move(pivot), num_pivots, num_per_leaf, page_size, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//        else if (program.get("--index") == "MVPTREE")
//        {
//            index = std::make_unique<gervLib::index::mvptree::MVPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, 2, 2, 4, 2, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//        else if (program.get("--index") == "OMNIKDTREE")
//        {
//            index = std::make_unique<gervLib::index::omni::OmniKdTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//        else if (program.get("--index") == "PMTREE")
//        {
//            index = std::make_unique<gervLib::index::pmtree::PMTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//        else if (program.get("--index") == "SPBTREE")
//        {
//            index = std::make_unique<gervLib::index::spbtree::SPBTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, num_bins, page_size, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//        else if (program.get("--index") == "LC")
//        {
//            index = std::make_unique<gervLib::index::lc::LC<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, path.string());
//        }
//
//        //index->saveIndex();
//
//    }
//    else if (program.get("operation") == "CALCULATE_LID")
//    {
//        std::cout << "Calculate LID" << std::endl;
//    }
//    else if (program.get("operation") == "QUERY")
//    {
//        std::cout << "Query" << std::endl;
//    }
//    else
//    {
//        throw std::invalid_argument("Invalid operation");
//    }
//
//
//    return 0;
//
//
//}



//#include "Configure.h"
//#include "Utils.h"
//#include "Dataset.h"
////#include "EditDistance.h"
//#include "EuclideanDistance.h"
//#include "RandomPivots.h"
//#include "KmedoidsPivots.h"
//#include "PageManager.h"
//#include "SequentialScan.h"
//#include <cassert>
//#include "LAESA.h"
//#include "VPTree.h"
//#include "MVPTree.h"
//#include "OmniKdTree.h"
////#include "KdTree.h"
//#include "PMTree.h"
//#include "SPBTree.h"
//#include "LC.h"
//#include "ISPivots.h"
//#include "FFTPivots.h"
//#include "DatasetUtils.h"
//
//using namespace gervLib::index;
//using namespace gervLib::configure;
//using namespace gervLib::utils;
//using namespace gervLib::dataset;
//using namespace gervLib::distance;
//using namespace gervLib::pivots;
//using namespace gervLib::memory;
//using namespace gervLib::query;
//
//
//int main(int argc, char **argv)
//{
//
//    gervLib::configure::configure();
//    std::cout << std::boolalpha;
//
//    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ","),
//            data2 = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ","),
//            test = std::make_unique<Dataset<size_t, double>>("../data/mnist_7k.csv", ",");
//
//    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
//            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
//
//    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
//    pvt->setSampleSize(1.0);
//    pvt->setSeed(449462);
//
//    std::unique_ptr<Index<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 12, 50, 0, false, false, true, true, "VP_SERI_TEST");
////    std::cout << *vp << std::endl;
////    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4096, 2, 2, 4, 2, false, true, true, true);
////    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 4096, false, true, true);
////    std::unique_ptr<pmtree::PMTree<size_t, double>> pm = std::make_unique<pmtree::PMTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true);
////    std::unique_ptr<spbtree::SPBTree<size_t, double>> spb = std::make_unique<spbtree::SPBTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 50, 4096, false, true, true);
////    std::unique_ptr<lc::LC<size_t, double>> lc = std::make_unique<lc::LC<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 0, true, true, true, true);
//
//    vp->saveIndex();
//
//    std::unique_ptr<Index<size_t, double>> vp2 = std::make_unique<vptree::VPTree<size_t, double>>("VP_SERI_TEST");
//
//    std::cout << vp->isEqual(vp2) << std::endl;
//
////    splitByLID(data2, dist2, vp, 100, "cities", "../data");
////    splitTrainTest<size_t, double>("../data/cities_norm.csv", ",", 0.8, 176);
//
////    std::vector<gervLib::query::ResultEntry<size_t>> res1 = pm->kNNIncremental(test->getElement(0), 100, true, true);
////    std::vector<gervLib::query::ResultEntry<size_t>> res2 = vp->kNNIncremental(test->getElement(0), 100, true, true);
////    std::cout << *lc << std::endl;
////
////    std::unique_ptr<u_char[]> serialized = lc->serialize();
////    std::unique_ptr<Index<size_t, double>> lc2 = std::make_unique<lc::LC<size_t, double>>();
////    lc2->deserialize(std::move(serialized));
////
////    std::cout << lc->isEqual(lc2) << std::endl;
//
////    std::unique_ptr<Index<size_t, double>> spb = std::make_unique<spbtree::SPBTree<size_t, double, mpz_class>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 5, 4096, false, false, true);
//
////    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
////    std::unique_ptr<spbtree::SPBTree<size_t, double>> spb2 = std::make_unique<spbtree::SPBTree<size_t, double>>();
////    std::unique_ptr<u_char[]> serialized = spb->serialize();
////    spb2->deserialize(std::move(serialized));
//
////    std::cout << *spb << "\n\n";
//
////    for(size_t i = 0; i < test->getCardinality(); i++)
////    {
////
////        std::vector<gervLib::query::ResultEntry<size_t>> res1 = lc2->kNNIncremental(test->getElement(i), 100, true, true);
////        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true, true);
////
////        for(size_t j = 0; j < res1.size(); j++)
////        {
////            if (res1[j].getDistance() != res2[j].getDistance()) {
////                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
////                //throw std::runtime_error("Error");
////            }
////
//////            std::cout << spb->getPrunning() << std::endl;
////
////        }
////    }
//
////    std::vector<gervLib::query::ResultEntry<size_t>> res1 = lc->kNNIncremental(test->getElement(0), 5, true, true);
////    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(0), 5, true);
////
////    for(size_t j = 0; j < res1.size(); j++)
////    {
////        std::cout << res2[j] << "\t" << res1[j] << std::endl;
////    }
//
////    size_t id = 2;
////    std::vector<gervLib::query::ResultEntry<size_t>> res1 = spb->kNNIncremental(test->getElement(id), 10, true);
////    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(id), 10, true);
////
////    for(size_t j = 0; j < res1.size(); j++)
////    {
////        std::cout << res1[j] << "\t" << res2[j] << std::endl;
////        if (res1[j].getDistance() != res2[j].getDistance()) {
////                //assert(res1[j].getDistance() == res2[j].getDistance());
////            std::cout << "Error: " << res1[j].getDistance() << " != " << res2[j].getDistance() << std::endl;
////            throw std::runtime_error("Error");
////        }
////    }
////
////    std::vector<gervLib::query::ResultEntry<size_t>> res3 = spb->kNNIncremental(test->getElement(id), 10, true);
////    std::vector<gervLib::query::ResultEntry<size_t>> res4 = sc->kNN(test->getElement(id), 10, true);
////
////    for(size_t j = 0; j < res3.size(); j++)
////    {
////        std::cout << res3[j] << "\t" << res4[j] << std::endl;
////        if (res3[j].getDistance() != res4[j].getDistance()) {
////            //assert(res1[j].getDistance() == res2[j].getDistance());
////            std::cout << "Error: " << res3[j].getDistance() << " != " << res4[j].getDistance() << std::endl;
////            throw std::runtime_error("Error");
////        }
////    }
//
//
////    std::unique_ptr<omni::OmniKdTree<size_t, double>> omni = std::make_unique<omni::OmniKdTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test11");
////    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
////
////    for(size_t i = 0; i < test->getCardinality(); i++)
////    {
////        std::vector<gervLib::query::ResultEntry<size_t>> res1 = omni->kNNIncremental(test->getElement(i), 100, true);
////        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
////
////        for(size_t j = 0; j < res1.size(); j++)
////        {
////            if (res1[j].getDistance() != res2[j].getDistance()) {
////                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
////                //throw std::runtime_error("Error");
////            }
//////            std::cout << omni->getPrunning() << std::endl;
////        }
////    }
//
////    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 5, 4096);
////    std::unique_ptr<u_char[]> serialized = mvp->serialize();
////    std::unique_ptr<Index<size_t, double>> mvp2 = std::make_unique<mvptree::MVPTree<size_t, double>>();
////    mvp2->deserialize(std::move(serialized));
////
////    std::cout << mvp->isEqual(mvp2) << std::endl;
//
////    std::unique_ptr<Dataset<size_t, double>> data1 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
////            data2 = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ","),
////            test = std::make_unique<Dataset<size_t, double>>("../data/cities_norm.csv", ",");
////    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> dist1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
////            dist2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
////
////    auto pvt = std::make_unique<RandomPivots<size_t, double>>();
////
////    std::unique_ptr<mvptree::MVPTree<size_t, double>> mvp = std::make_unique<mvptree::MVPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2,
////            500, 4096, 2, 2, 4, 2, false, true, true, false, "tmp_unit_test11");
////    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
////
////    for(size_t i = 0; i < test->getCardinality(); i++)
////    {
////        std::vector<gervLib::query::ResultEntry<size_t>> res1 = mvp->kNNIncremental(test->getElement(i), 100, true);
////        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
////
////        for(size_t j = 0; j < res1.size(); j++)
////            assert(res1[j].getDistance() == res2[j].getDistance());
////
////        std::cout << mvp->getPrunning() << std::endl;
////
////    }
////
////    gervLib::utils::deleteDirectory("tmp_unit_test11");
////    gervLib::utils::deleteDirectory("tmp_unit_test12");
//
//
////    std::unique_ptr<vptree::VPTree<size_t, double>> vp = std::make_unique<vptree::VPTree<size_t, double>>(std::move(data1), std::move(dist1), std::move(pvt), 2, 50, 8000, false, true, true, "tmp_unit_test11");
////    std::unique_ptr<SequentialScan<size_t, double>> sc = std::make_unique<SequentialScan<size_t, double>>(std::move(data2), std::move(dist2), "tmp_unit_test12");
////
////    for(size_t i = 0; i < test->getCardinality(); i++)
////    {
////        std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(i), 100, true);
////        std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(i), 100, true);
////
////        for(size_t j = 0; j < res1.size(); j++)
////        {
////            if (res1[j].getDistance() != res2[j].getDistance()) {
////                std::cout << "Error index " << i << ": " << res1[j].getDistance() << " != " << res2[j].getDistance()<< std::endl;
////                throw std::runtime_error("Error");
////            }
////            std::cout << vp->getPrunning() << std::endl;
////        }
////    }
//
//////    std::vector<gervLib::query::ResultEntry<size_t>> res1 = vp->kNNIncremental(test->getElement(47), 5, true);
//////    std::vector<gervLib::query::ResultEntry<size_t>> res2 = sc->kNN(test->getElement(6), 100, true);
////
//////    for(size_t j = 0; j < res1.size(); j++)
//////    {
//////        if (res1[j].getDistance() != res2[j].getDistance()) {
//////                //assert(res1[j].getDistance() == res2[j].getDistance());
//////            std::cout << "Error index " << 0 << ": " << res1[j].getDistance() << " != " << res2[j].getDistance() << std::endl;
//////            throw std::runtime_error("Error");
//////        }
//////    }
////
////    gervLib::utils::deleteDirectory("tmp_unit_test11");
////    gervLib::utils::deleteDirectory("tmp_unit_test12");
////
//    return 0;
//
//
//}
//
#include <iostream>
#include "Indexes.h"
#include "Pivots.h"
#include "Query.h"

/* KeyWords:
 * -INDEX:
 * -DATASET_NAME:
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
 * -NUM_PER_LEAF: Default | DATASET_TRAIN | * 0.01
 * -NUM_BINS: Default 256
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
        throw std::invalid_argument("Invalid number of arguments");
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

        if (args.find("-DATASET_NAME") == args.end() || args.find("-INDEX") == args.end() || args.find("-DATASET_TRAIN") == args.end() || args.find("-DATASET_TEST") == args.end() || args.find("-PIVOT_TYPE") == args.end())
        {

            throw std::invalid_argument("Index, dataset name, dataset train, dataset test or pivot type not specified");

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

            num_bins = 256;

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

            use_laesa = false;

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

        std::string prefix = args["-INDEX"] + "_" + args["-DATASET_NAME"] + "_" + args["-PIVOT_TYPE"];
        std::string indexPath = gervLib::utils::generatePathByPrefix(gervLib::configure::baseOutputPath, prefix, true);

        std::cout << "Building index...\n";

        //Index
        if (args["-INDEX"] == "VPTREE")
        {
            index = std::make_unique<gervLib::index::vptree::VPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }
        else if (args["-INDEX"] == "MVPTREE")
        {
            index = std::make_unique<gervLib::index::mvptree::MVPTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, 2, 2, 4, 2, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }
        else if (args["-INDEX"] == "OMNIKDTREE")
        {
            index = std::make_unique<gervLib::index::omni::OmniKdTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }
        else if (args["-INDEX"] == "PMTREE")
        {
            index = std::make_unique<gervLib::index::pmtree::PMTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }
        else if (args["-INDEX"] == "SPBTREE")
        {
            index = std::make_unique<gervLib::index::spbtree::SPBTree<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, num_bins, page_size, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }
        else if (args["-INDEX"] == "LC")
        {
            index = std::make_unique<gervLib::index::lc::LC<size_t, double>>(std::move(dataset_train), std::move(distance_function), std::move(pivot), num_pivots, num_per_leaf, page_size, store_pivot_leaf, store_directory_node, store_leaf_node, use_laesa, indexPath);
        }

//        //Queries
//        size_t max_queries = std::min(dataset_test->getCardinality(), num_query);
//        std::vector<ResultEntry<size_t>> results;
//
//        for (size_t i = 0; i < dataset_test->getCardinality(); i++)
//        {
//            index->kNNIncremental(dataset_test->getElement(i), k_max, false, true, results);
//
//            if ((i + 1) % num_queries_per_file == 0)
//                index->generateIndexFiles(false, true);
//        }

    }

    //PCA, KMEDOIDS, MAXVARIANCE

    return 0;

}