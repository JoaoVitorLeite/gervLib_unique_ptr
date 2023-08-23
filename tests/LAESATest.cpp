//
// Created by joaovictor on 06/08/23.
//

#include "LAESA.h"
#include "SequentialScan.h"

using namespace gervLib::index;
using namespace gervLib::dataset;
using namespace gervLib::distance;
using namespace gervLib::query;
using namespace gervLib::pivots;

int test1()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " "),
                                             data1 = std::make_unique<Dataset<size_t, double>>(*data), data2 = std::make_unique<Dataset<size_t, double>>(*data);

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
                                                                        df2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data1), std::move(df1), "tmp_unit_test4");

    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<KmedoidsPivots<size_t, double>>();
    std::unique_ptr<Index<size_t, double>> laesa = std::make_unique<LAESA<size_t, double>>(std::move(data2), std::move(df2), std::move(pivots), 2, "tmp_unit_test5");

    std::vector<ResultEntry<size_t>> resultsLaesa, resultsSequentialScan;

    for(size_t i = 0; i < data->getCardinality(); i++)
    {

        resultsLaesa = laesa->kNN(data->getElement(i), 10, true);
        resultsSequentialScan = sequentialScan->kNN(data->getElement(i), 10, true);

        for(size_t j = 0; j < resultsLaesa.size(); j++)
        {

            if(resultsLaesa[j] != resultsSequentialScan[j])
            {

                std::cout << "Error: LAESA and SequentialScan results differ for element " << i << std::endl;
                return 1;

            }

        }

        resultsLaesa.clear();
        resultsSequentialScan.clear();

    }

    gervLib::utils::deleteDirectory("tmp_unit_test4");
    gervLib::utils::deleteDirectory("tmp_unit_test5");

    return 0;

}

int test2()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " "),
            data1 = std::make_unique<Dataset<size_t, double>>(*data), data2 = std::make_unique<Dataset<size_t, double>>(*data);

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            df2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data1), std::move(df1), "tmp_unit_test6");

    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<KmedoidsPivots<size_t, double>>();
    std::unique_ptr<Index<size_t, double>> laesa = std::make_unique<LAESA<size_t, double>>(std::move(data2), std::move(df2), std::move(pivots), 2, "tmp_unit_test7");

    std::vector<ResultEntry<size_t>> resultsLaesa, resultsSequentialScan;

    for(size_t i = 0; i < data->getCardinality(); i++)
    {

        resultsLaesa = laesa->kNN(data->getElement(i), 10, true);
        resultsSequentialScan = sequentialScan->kNN(data->getElement(i), 10, true);

        for(size_t j = 0; j < resultsLaesa.size(); j++)
        {

            if(resultsLaesa[j] != resultsSequentialScan[j])
            {

                std::cout << "Error: LAESA and SequentialScan results differ for element " << i << std::endl;
                return 1;

            }

        }

        resultsLaesa.clear();
        resultsSequentialScan.clear();

    }

    gervLib::utils::deleteDirectory("tmp_unit_test6");
    gervLib::utils::deleteDirectory("tmp_unit_test7");

    return 0;

}

int test3()
{


    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/Dataset1.csv", " ");
    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> distanceFunction = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();
    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<KmedoidsPivots<size_t, double>>();
    std::unique_ptr<Index<size_t, double>> index = std::make_unique<LAESA<size_t, double>>(std::move(data), std::move(distanceFunction), std::move(pivots), 2, "tmp_unit_test8");

    index->saveIndex();

    std::unique_ptr<Index<size_t, double>> index2 = std::make_unique<LAESA<size_t, double>>("tmp_unit_test8", "");

    assert(index->isEqual(index2));

    gervLib::utils::deleteDirectory("tmp_unit_test8");

    return 0;

}

int test4()
{

    std::unique_ptr<Dataset<size_t, double>> data = std::make_unique<Dataset<size_t, double>>("../../data/cities_norm.csv", ","),
            data1 = std::make_unique<Dataset<size_t, double>>(*data), data2 = std::make_unique<Dataset<size_t, double>>(*data);

    std::unique_ptr<DistanceFunction<BasicArrayObject<size_t, double>>> df1 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>(),
            df2 = std::make_unique<EuclideanDistance<BasicArrayObject<size_t, double>>>();

    std::unique_ptr<SequentialScan<size_t, double>> sequentialScan = std::make_unique<SequentialScan<size_t, double>>(std::move(data1), std::move(df1), "tmp_unit_test4");

    std::unique_ptr<Pivot<size_t, double>> pivots = std::make_unique<KmedoidsPivots<size_t, double>>();
    std::unique_ptr<Index<size_t, double>> laesa = std::make_unique<LAESA<size_t, double>>(std::move(data2), std::move(df2), std::move(pivots), 2, "tmp_unit_test5");

    std::vector<ResultEntry<size_t>> resultsLaesa, resultsSequentialScan;

    for(size_t i = 0; i < data->getCardinality(); i++)
    {

        resultsLaesa = laesa->kNN(data->getElement(i), 100, true);
        resultsSequentialScan = sequentialScan->kNN(data->getElement(i), 100, true);

        for(size_t j = 0; j < resultsLaesa.size(); j++)
        {

            if(resultsLaesa[j] != resultsSequentialScan[j])
            {

                std::cout << "Error: LAESA and SequentialScan results differ for element " << i << std::endl;
                return 1;

            }

        }

        resultsLaesa.clear();
        resultsSequentialScan.clear();

    }

    gervLib::utils::deleteDirectory("tmp_unit_test4");
    gervLib::utils::deleteDirectory("tmp_unit_test5");

    return 0;

}

int main(int argc, char *argv[])
{

    gervLib::configure::configure();
    std::cout << std::boolalpha;

    int res = 0;

    res += test1();
    res += test2();
    res += test3();
    res += test4();

    return res == 0 ? 0 : 1;

}