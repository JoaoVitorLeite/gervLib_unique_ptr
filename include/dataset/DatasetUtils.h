//
// Created by joaoleite on 9/12/23.
//

#ifndef GERVLIB_DATASETUTILS_H
#define GERVLIB_DATASETUTILS_H

#include "Dataset.h"
#include "DistanceFunction.h"
#include "Index.h"
#include "Utils.h"

namespace gervLib::dataset
{

    template <typename O, typename T>
    void splitByLID(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                    std::unique_ptr<gervLib::distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df,
                    std::unique_ptr<gervLib::index::Index<O, T>>& index, size_t k, std::string fileName, std::string outputPath)
    {

        std::vector<query::ResultEntry<O>> result;
        double lastDist, lid;
        std::vector<std::pair<double, size_t>> lids;

        for (size_t i = 0; i < dataset->getCardinality(); i++)
        {
            result = index->kNNIncremental(dataset->getElement(i), k+1, false, false);
            result.erase(result.begin());

            lastDist = result.back().getDistance();
            lid = 0.0;

            for (size_t j = 0; j < k; j++)
            {
                lid += log(result[j].getDistance() / lastDist);
            }

            lid /= static_cast<double>(k);
            lid = -1.0/lid;

            lids.emplace_back(lid, i);
            result.clear();
        }

        std::sort(lids.begin(), lids.end());

        std::filesystem::path path(outputPath);
        path /= fileName;
        std::string q1FileName = path.string() + "_Q1.csv",
                    q2FileName = path.string() + "_Q2.csv",
                    q3FileName = path.string() + "_Q3.csv",
                    q4FileName = path.string() + "_Q4.csv";

        std::ofstream q1File(q1FileName), q2File(q2FileName), q3File(q3FileName), q4File(q4FileName);
        size_t n = lids.size();

        size_t lowerHalf = n / 2;
        double Q1 = (n % 2 == 0) ? (lids[lowerHalf/2 - 1].first + lids[lowerHalf/2].first) / 2.0 : lids[lowerHalf/2].first,
               Q2 = (n % 2 == 0) ? (lids[n/2 - 1].first + lids[n/2].first) / 2.0 : lids[n/2].first,
               Q3 = (n % 2 == 0) ? (lids[n/2 + lowerHalf/2 - 1].first + lids[n/2 + lowerHalf/2].first) / 2.0 : lids[n/2 + lowerHalf/2].first;

        for (size_t i = 0; i < n; i++)
        {

            if (lids[i].first <= Q1)
            {
                for (size_t x = 0; x < dataset->getElement(lids[i].second).size(); x++)
                {
                    if (x != dataset->getElement(lids[i].second).size() - 1)
                        q1File << dataset->getElement(lids[i].second).operator[](x) << ",";
                    else
                        q1File << dataset->getElement(lids[i].second).operator[](x) << std::endl;
                }
            }
            else if (lids[i].first <= Q2)
            {
                for (size_t x = 0; x < dataset->getElement(lids[i].second).size(); x++)
                {
                    if (x != dataset->getElement(lids[i].second).size() - 1)
                        q2File << dataset->getElement(lids[i].second).operator[](x) << ",";
                    else
                        q2File << dataset->getElement(lids[i].second).operator[](x) << std::endl;
                }
            }
            else if (lids[i].first <= Q3)
            {
                for (size_t x = 0; x < dataset->getElement(lids[i].second).size(); x++)
                {
                    if (x != dataset->getElement(lids[i].second).size() - 1)
                        q3File << dataset->getElement(lids[i].second).operator[](x) << ",";
                    else
                        q3File << dataset->getElement(lids[i].second).operator[](x) << std::endl;
                }
            }
            else
            {
                for (size_t x = 0; x < dataset->getElement(lids[i].second).size(); x++)
                {
                    if (x != dataset->getElement(lids[i].second).size() - 1)
                        q4File << dataset->getElement(lids[i].second).operator[](x) << ",";
                    else
                        q4File << dataset->getElement(lids[i].second).operator[](x) << std::endl;
                }
            }

        }

        lids.clear();

        q1File.close();
        q2File.close();
        q3File.close();
        q4File.close();

    }

    template <typename O, typename T>
    void calculateLID(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                    std::unique_ptr<gervLib::distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df,
                    std::unique_ptr<gervLib::index::Index<O, T>>& index, size_t k, std::string fileName, std::string outputPath)
    {

        std::vector<query::ResultEntry<O>> result;
        double lastDist, lid;
        std::vector<std::pair<double, size_t>> lids;

        for (size_t i = 0; i < dataset->getCardinality(); i++)
        {
            result = index->kNNIncremental(dataset->getElement(i), k+1, false, false);
            result.erase(result.begin());

            lastDist = result.back().getDistance();
            lid = 0.0;

            for (size_t j = 0; j < k; j++)
            {
                lid += log(result[j].getDistance() / lastDist);
            }

            lid /= static_cast<double>(k);
            lid = -1.0/lid;

            lids.emplace_back(lid, dataset->getElement(i).getOID());
            result.clear();
        }

        std::filesystem::path path(outputPath);
        path /= fileName;
        std::string LIDFileName = path.string() + "_lid.csv";

        std::ofstream LIDFile(LIDFileName);
        size_t n = lids.size();

        for (size_t i = 0; i < n; i++)
        {

            LIDFile << lids[i].first << "," << lids[i].second << std::endl;

        }

        lids.clear();

        LIDFile.close();

    }

    template <typename O, typename T>
    void splitTrainTest(std::string datasetPath, std::string separator, std::variant<size_t, double> _trainSize, size_t seed)
    {

        size_t trainSize, testSize;
        std::unique_ptr<Dataset<O, T>> dataset = std::make_unique<Dataset<O, T>>(datasetPath, separator);

        if (std::holds_alternative<size_t>(_trainSize))
        {
            trainSize = std::get<size_t>(_trainSize);
            testSize = dataset->getCardinality() - trainSize;
        }
        else
        {
            double trainSizePerc = std::get<double>(_trainSize);
            trainSize = static_cast<size_t>(std::ceil(dataset->getCardinality() * trainSizePerc));
            testSize = dataset->getCardinality() - trainSize;
        }

        utils::Random<size_t> random(seed);
        std::vector<size_t> indexes(dataset->getCardinality());
        std::iota(indexes.begin(), indexes.end(), 0);
        std::shuffle(indexes.begin(), indexes.end(), random.generator);

        std::filesystem::path path(datasetPath);
        std::string fileName = path.stem().string(),
                    trainFileName = path.replace_filename(fileName + "_train.csv").string(),
                    testFileName = path.replace_filename(fileName + "_test.csv").string();

        std::ofstream trainFile(trainFileName), testFile(testFileName);

        for (size_t i = 0; i < trainSize; i++)
        {
            for (size_t j = 0; j < dataset->getElement(indexes[i]).size(); j++)
            {
                if (j != dataset->getElement(indexes[i]).size() - 1)
                    trainFile << dataset->getElement(indexes[i]).operator[](j) << ",";
                else
                    trainFile << dataset->getElement(indexes[i]).operator[](j) << std::endl;
            }
        }

        for (size_t i = trainSize; i < dataset->getCardinality(); i++)
        {
            for (size_t j = 0; j < dataset->getElement(indexes[i]).size(); j++)
            {
                if (j != dataset->getElement(indexes[i]).size() - 1)
                    testFile << dataset->getElement(indexes[i]).operator[](j) << ",";
                else
                    testFile << dataset->getElement(indexes[i]).operator[](j) << std::endl;
            }
        }

        trainFile.close();
        testFile.close();

    }

}

#endif //GERVLIB_DATASETUTILS_H
