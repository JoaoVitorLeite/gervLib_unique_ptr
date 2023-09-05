//
// Created by joaovictor on 17/08/23.
//

#ifndef GERVLIB_INDEX_H
#define GERVLIB_INDEX_H

#include "Dataset.h"
#include "DistanceFactory.h"
#include "PivotFactory.h"
#include "PageManager.h"
#include "Query.h"

namespace gervLib::index
{

    enum INDEX_TYPE {SEQUENTIAL_SCAN_t, LAESA_t, VPTREE_t, MVPTREE, KDTREE, OMNIKDTREE, PMTREE, SPBTREE, UNKNOWN};

    enum MEMORY_STATUS {IN_MEMORY, IN_DISK, NONE};

    std::map<INDEX_TYPE, std::string> indexTypeMap = {
            {SEQUENTIAL_SCAN_t, "SEQUENTIAL_SCAN"},
            {LAESA_t, "LAESA"},
            {VPTREE_t, "VPTREE"},
            {MVPTREE, "MVPTREE"},
            {KDTREE, "KDTREE"},
            {OMNIKDTREE, "OMNIKDTREE"},
            {PMTREE, "PMTREE"},
            {SPBTREE, "SPBTREE"}
    };

    std::map<std::string, INDEX_TYPE> indexTypeMapReverse = {
            {"SEQUENTIAL_SCAN", SEQUENTIAL_SCAN_t},
            {"LAESA", LAESA_t},
            {"VPTREE", VPTREE_t},
            {"MVPTREE", MVPTREE},
            {"KDTREE", KDTREE},
            {"OMNIKDTREE", OMNIKDTREE},
            {"PMTREE", PMTREE},
            {"SPBTREE", SPBTREE}
    };

    std::map<MEMORY_STATUS, std::string> memoryStatusMap = {
            {IN_MEMORY, "IN_MEMORY"},
            {IN_DISK, "IN_DISK"},
            {NONE, "NONE"}
    };

    std::map<std::string, MEMORY_STATUS> memoryStatusMapReverse = {
            {"IN_MEMORY", IN_MEMORY},
            {"IN_DISK", IN_DISK},
            {"NONE", NONE}
    };

    template <typename O, typename T>
    class Index : public gervLib::serialize::Serialize
    {

    protected:
        std::unique_ptr<dataset::Dataset<O, T>> dataset;
        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O,T>>> distanceFunction;
        std::unique_ptr<pivots::Pivot<O,T>> pivots;
        std::unique_ptr<memory::PageManager<O>> pageManager;
        size_t pageSize{}, prunning{}, leafNodeAccess{};
        INDEX_TYPE indexType = UNKNOWN;
        std::string indexName, indexFolder, buildFile, experimentFile;

    protected:
        virtual std::string headerBuildFile() = 0;

        virtual std::string headerExperimentFile() = 0;

        void saveStatistics(std::initializer_list<std::string> list)
        {

            std::ofstream file(experimentFile, std::ios::app);

            if(file.is_open())
            {

                std::string accum = std::accumulate(list.begin(), list.end(), std::string{}, [](const std::string& a, const std::string& b) {
                    return a.empty() ? b : a + "," + b;
                });

                file << accum << std::endl;

            }
            else
            {
                throw std::runtime_error("Index::saveStatistics(): Error opening file");
            }

            file.close();

        }

        void saveResultToFile(std::vector<gervLib::query::ResultEntry<O>> &result, dataset::BasicArrayObject<O, T>& query, std::string queryName, std::string expt_id, std::initializer_list<std::string> statsList)
        {

            size_t total = sizeof(size_t)*3 + query.getSerializedSize(), offset = 0, sz = result.size();

            for(auto &entry : result)
                total += sizeof(size_t) + entry.getSerializedSize();

            std::unique_ptr<u_char[]> data(new u_char[total]);

            std::memcpy(data.get() + offset, &total, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            for(auto &entry : result) {

                sz = entry.getSerializedSize();
                std::memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> entryData = entry.serialize();
                std::memcpy(data.get() + offset, entryData.get(), sz);
                offset += sz;
                entryData.reset();

            }

            sz = query.getSerializedSize();
            std::memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> queryData = query.serialize();
            std::memcpy(data.get() + offset, queryData.get(), sz);
            offset += sz;
            queryData.reset();

            std::string fileName = this->indexName + "_" + queryName + "_" + expt_id + ".query";
            std::filesystem::path path(this->indexFolder);
            path /= fileName;

            std::ofstream file(path.string(), std::ios::binary);

            if (file.is_open()) {
                file.write((char *) data.get(), total);
                file.close();
            }
            else
                throw std::runtime_error("Index::saveResultToFile(): Error opening file");

            saveStatistics(statsList);

            data.reset();

        }

    public:
        Index() = default;

        virtual ~Index()
        {
            if (dataset != nullptr)
                dataset.reset();

            if (distanceFunction != nullptr)
                distanceFunction.reset();

            if (pivots != nullptr)
                pivots.reset();

            if (pageManager != nullptr)
                pageManager.reset();
        }

        void generateIndexFiles(bool build, bool experiment)
        {

            if (indexFolder.empty())
                indexFolder = gervLib::utils::generatePathByPrefix(configure::baseOutputPath, indexName);

            utils::createFolderIfNotExists(indexFolder);

            if (build)
            {

                buildFile = utils::generateFileByPrefix(indexFolder, "build", true, ".txt");
                std::ofstream buildFileStream(buildFile);
                buildFileStream << headerBuildFile() << std::endl;
                buildFileStream.close();

            }

            if (experiment)
            {

                experimentFile = utils::generateFileByPrefix(indexFolder, "experiment", true, ".txt");
                std::ofstream experimentFileStream(experimentFile);
                experimentFileStream << headerExperimentFile() << std::endl;
                experimentFileStream.close();

            }

        }

        void saveIndex()
        {

            std::string path = gervLib::utils::generateFileByPrefix(this->indexFolder, this->indexName, true, ".index");
            std::ofstream file(path, std::ios::binary);

            if (file.is_open())
            {

                std::unique_ptr<u_char[]> data = serialize();
                size_t size = getSerializedSize();
                std::unique_ptr<u_char[]> wrt = std::make_unique<u_char[]>(sizeof(size_t) + size);

                std::memcpy(wrt.get(), &size, sizeof(size_t));
                std::memcpy(wrt.get() + sizeof(size_t), data.get(), size);

                file.write((char *)wrt.get(), sizeof(size_t) + size);

                data.reset();
                wrt.reset();

            }
            else
            {
                throw std::runtime_error("Index::saveIndex(): Error opening file");
            }

            file.close();

        }

        void loadIndex(std::string indexPath)
        {

            std::ifstream file(indexPath, std::ios::binary);

            if(file.is_open())
            {

                size_t size;
                std::unique_ptr<u_char[]> sizeChar = std::make_unique<u_char[]>(sizeof(size_t));
                file.read((char *)sizeChar.get(), sizeof(size_t));
                std::memcpy(&size, sizeChar.get(), sizeof(size_t));
                sizeChar.reset();

                std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(size);
                file.read((char *)data.get(), size);
                deserialize(std::move(data));

            }
            else
            {
                throw std::runtime_error("Index::loadIndex(): Error opening file");
            }

            file.close();

        }

        void loadIndex()
        {

            std::string indexPath = utils::getFileByPrefix(this->indexFolder, ".index");
            std::filesystem::path path(this->indexFolder);
            path /= indexPath;
            std::ifstream file(path.string(), std::ios::binary);

            if(file.is_open())
            {

                size_t size;
                std::unique_ptr<u_char[]> sizeChar = std::make_unique<u_char[]>(sizeof(size_t));
                file.read((char *)sizeChar.get(), sizeof(size_t));
                std::memcpy(&size, sizeChar.get(), sizeof(size_t));
                sizeChar.reset();

                std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(size);
                file.read((char *)data.get(), size);
                deserialize(std::move(data));

            }
            else
            {
                throw std::runtime_error("Index::loadIndex(): Error opening file");
            }

            file.close();

        }

        std::pair<gervLib::dataset::BasicArrayObject<O, T>, std::vector<gervLib::query::ResultEntry<O>>> loadResultFromFile(std::string fileName)
        {

            std::vector<gervLib::query::ResultEntry<O>> result;

            std::unique_ptr<u_char[]> dataSize = std::make_unique<u_char[]>(sizeof(size_t));
            size_t total, offset = 0, sz, num_elements;

            std::ifstream file(fileName, std::ios::binary);
            file.read((char *)dataSize.get(), sizeof(size_t));

            std::memcpy(&total, dataSize.get(), sizeof(size_t));
            dataSize.reset();

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(total);
            file.read((char *)data.get(), total);
            file.close();

            std::memcpy(&num_elements, data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            for(size_t i = 0; i < num_elements; i++)
            {

                std::memcpy(&sz, data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> entryData = std::make_unique<u_char[]>(sz);
                std::memcpy(entryData.get(), data.get() + offset, sz);
                offset += sz;

                gervLib::query::ResultEntry<O> entry;
                entry.deserialize(std::move(entryData));
                result.push_back(entry);

            }

            std::memcpy(&sz, data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> queryData = std::make_unique<u_char[]>(sz);
            std::memcpy(queryData.get(), data.get() + offset, sz);
            offset += sz;

            gervLib::dataset::BasicArrayObject<O, T> query;
            query.deserialize(std::move(queryData));

            data.reset();

            return std::make_pair(query, result);

        }

        [[nodiscard]] size_t getDistanceCount() const
        {
            return distanceFunction->getDistanceCount();
        }

        [[nodiscard]] INDEX_TYPE getIndexType() const
        {
            return indexType;
        }

        void setIndexType(INDEX_TYPE _indexType)
        {
            this->indexType = _indexType;
        }

        [[nodiscard]] size_t getPageSize() const
        {
            return pageSize;
        }

        [[nodiscard]] const std::string &getIndexName() const
        {
            return indexName;
        }

        void setIndexName(const std::string &_indexName)
        {
            indexName = _indexName;
        }

        [[nodiscard]] const std::string &getIndexFolder() const
        {
            return indexFolder;
        }

        void setIndexFolder(const std::string &_indexFolder)
        {
            indexFolder = _indexFolder;
        }

        [[nodiscard]] const std::string &getBuildFile() const
        {
            return buildFile;
        }

        [[nodiscard]] const std::string &getExperimentFile() const
        {
            return experimentFile;
        }

        [[nodiscard]] size_t getPrunning() const
        {
            return prunning;
        }

        [[nodiscard]] size_t getLeafNodeAccess() const
        {
            return leafNodeAccess;
        }

        std::unique_ptr<dataset::Dataset<O, T>>& getDataset()
        {
            return dataset;
        }

        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& getDistanceFunction()
        {
            return distanceFunction;
        }

        std::unique_ptr<pivots::Pivot<O, T>>& getPivots()
        {
            return pivots;
        }

        std::unique_ptr<memory::PageManager<O>>& getPageManager()
        {
            return pageManager;
        }

        void setDataset(std::unique_ptr<dataset::Dataset<O, T>> _dataset)
        {
            if (dataset != nullptr)
            {
                dataset->clear();
                dataset.reset();
            }
            dataset = std::move(_dataset);
        }

        void setDistanceFunction(std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df)
        {
            if (distanceFunction != nullptr)
            {
                distanceFunction.reset();
            }
            distanceFunction = std::move(df);
        }

        void setPivots(std::unique_ptr<pivots::Pivot<O, T>> pvt)
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }
            pivots = std::move(pvt);
        }

        void kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults, std::vector<gervLib::query::ResultEntry<O>>& results)
        {

            results.clear();
            results = this->kNNIncremental(query, k, saveResults);

        }

        void kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults, std::vector<gervLib::query::ResultEntry<O>>& results)
        {

            results.clear();
            results = this->kNN(query, k, saveResults);

        }

        bool operator==(std::unique_ptr<Index<O, T>>& other)
        {

            return isEqual(other);

        }

        bool operator!=(std::unique_ptr<Index<O, T>>& other)
        {

            return !isEqual(other);

        }

        // Virtual methods

        virtual void clear()
        {
            if (dataset != nullptr)
                dataset.reset();

            if (distanceFunction != nullptr)
                distanceFunction.reset();

            if (pivots != nullptr)
                pivots.reset();

            if (pageManager != nullptr)
                pageManager.reset();
        }

        virtual void print(std::ostream& os) const
        {

            os << "Index Name: " << indexName << std::endl;
            os << "Index Type: " << indexTypeMap.at(indexType) << std::endl;
            os << "Index Folder: " << indexFolder << std::endl;
            os << "Build File: " << buildFile << std::endl;
            os << "Experiment File: " << experimentFile << std::endl;
            os << "Page Size: " << pageSize << std::endl;

            if (dataset != nullptr)
                os << "Dataset: " << std::endl << *dataset << std::endl;
            else
                os << "Dataset: " << std::endl << "NULL" << std::endl;

            if (distanceFunction != nullptr)
                os << "Distance Function: " << std::endl << *distanceFunction << std::endl;
            else
                os << "Distance Function: " << std::endl << "NULL" << std::endl;

            if (pivots != nullptr)
                os << "Pivots: " << std::endl << *pivots << std::endl;
            else
                os << "Pivots: " << std::endl << "NULL" << std::endl;

            if (pageManager != nullptr)
                os << "Page Manager: " << std::endl << *pageManager << std::endl;
            else
                os << "Page Manager: " << std::endl << "NULL" << std::endl;

        }

        virtual bool isEqual(std::unique_ptr<Index<O, T>>& other)
        {

            if (other == nullptr)
                return false;

            if (indexType != other->getIndexType())
                return false;

            if ((dataset == nullptr && other->getDataset() != nullptr) || (dataset != nullptr && other->getDataset() == nullptr))
                return false;

            if (dataset != nullptr && other->getDataset() != nullptr && !dataset->isEqual(*other->getDataset()))
                return false;

            if ((distanceFunction == nullptr && other->getDistanceFunction() != nullptr) || (distanceFunction != nullptr && other->getDistanceFunction() == nullptr))
                return false;

            if (distanceFunction != nullptr && other->getDistanceFunction() != nullptr && !distanceFunction->isEqual(other->getDistanceFunction()))
                return false;

            if ((pivots == nullptr && other->getPivots() != nullptr) || (pivots != nullptr && other->getPivots() == nullptr))
                return false;

            if (pivots != nullptr && other->getPivots() != nullptr && !pivots->isEqual(other->getPivots()))
                return false;

            return true;

        }

        virtual void buildIndex() = 0;

        virtual std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults)
        {
            throw std::runtime_error("Index::kNN(): Method not implemented");
        }

        virtual std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults)
        {
            throw std::runtime_error("Index::kNNIncremental(): Method not implemented");
        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;
            std::string str;

            std::memcpy(data.get() + offset, &pageSize, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, &prunning, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, &leafNodeAccess, sizeof(size_t));
            offset += sizeof(size_t);

            size = indexName.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, indexName.c_str(), size);
            offset += size;

            size = indexFolder.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, indexFolder.c_str(), size);
            offset += size;

            size = buildFile.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, buildFile.c_str(), size);
            offset += size;

            size = experimentFile.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, experimentFile.c_str(), size);
            offset += size;

            str = indexTypeMap[indexType];
            size = str.size();

            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, str.c_str(), size);
            offset += size;

            if (dataset != nullptr) {
                std::unique_ptr<u_char[]> datasetData = dataset->serialize();
                size = dataset->getSerializedSize();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
                std::memcpy(data.get() + offset, datasetData.get(), size);
                offset += size;
                datasetData.reset();
            }
            else
            {
                size = 0;
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (distanceFunction != nullptr) {

                size = distanceFunction->getSerializedSize();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                str = distance::distanceTypeMap[distanceFunction->getDistanceType()];
                size_t str_len = str.size();

                std::memcpy(data.get() + offset, &str_len, sizeof(size_t));
                offset += sizeof(size_t);

                std::memcpy(data.get() + offset, str.c_str(), str_len);
                offset += str_len;

                std::unique_ptr<u_char[]> distanceFunctionData = distanceFunction->serialize();
                std::memcpy(data.get() + offset, distanceFunctionData.get(), size);
                offset += size;

                distanceFunctionData.reset();

            }
            else
            {
                size = 0;
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (pivots != nullptr) {

                size = pivots->getSerializedSize();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                str = pivots::PIVOT_TYPE2STR[pivots->getPivotType()];
                size_t str_len = str.size();

                std::memcpy(data.get() + offset, &str_len, sizeof(size_t));
                offset += sizeof(size_t);

                std::memcpy(data.get() + offset, str.c_str(), str_len);
                offset += str_len;

                std::unique_ptr<u_char[]> pivotsData = pivots->serialize();
                std::memcpy(data.get() + offset, pivotsData.get(), size);
                offset += size;

                pivotsData.reset();

            }
            else
            {
                size = 0;
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (pageManager != nullptr) {

                size = pageManager->getSerializedSize();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> pageManagerData = pageManager->serialize();
                std::memcpy(data.get() + offset, pageManagerData.get(), size);
                offset += size;

                pageManagerData.reset();

            }
            else
            {
                size = 0;
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, size;
            std::string str;

            std::memcpy(&pageSize, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(&prunning, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(&leafNodeAccess, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            indexName.resize(size);
            std::memcpy(&indexName[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            indexFolder.resize(size);
            std::memcpy(&indexFolder[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            buildFile.resize(size);
            std::memcpy(&buildFile[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            experimentFile.resize(size);
            std::memcpy(&experimentFile[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            str.resize(size);
            std::memcpy(&str[0], _data.get() + offset, size);
            offset += size;

            indexType = indexTypeMapReverse[str];

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size != 0) {

                std::unique_ptr<u_char[]> datasetData = std::make_unique<u_char[]>(size);
                std::memcpy(datasetData.get(), _data.get() + offset, size);
                offset += size;

                if (dataset != nullptr)
                    dataset.reset();

                dataset = std::make_unique<dataset::Dataset<O, T>>();
                dataset->deserialize(std::move(datasetData));

            }

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size != 0) {

                size_t str_len;
                std::memcpy(&str_len, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                str.resize(str_len);
                std::memcpy(&str[0], _data.get() + offset, str_len);
                offset += str_len;

                if (distanceFunction != nullptr)
                    distanceFunction.reset();

                distanceFunction = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(distance::distanceNameMap[str]);

                std::unique_ptr<u_char[]> distanceFunctionData = std::make_unique<u_char[]>(size);
                std::memcpy(distanceFunctionData.get(), _data.get() + offset, size);
                offset += size;

                distanceFunction->deserialize(std::move(distanceFunctionData));

            }

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size != 0) {

                size_t str_len;
                std::memcpy(&str_len, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                str.resize(str_len);
                std::memcpy(&str[0], _data.get() + offset, str_len);
                offset += str_len;

                if (pivots != nullptr)
                    pivots.reset();

                pivots = pivots::PivotFactory<O, T>::createPivot(pivots::STR2PIVOT_TYPE[str]);

                std::unique_ptr<u_char[]> pivotsData = std::make_unique<u_char[]>(size);
                std::memcpy(pivotsData.get(), _data.get() + offset, size);
                offset += size;

                pivots->deserialize(std::move(pivotsData));

            }

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size != 0) {

                std::unique_ptr<u_char[]> pageManagerData = std::make_unique<u_char[]>(size);
                std::memcpy(pageManagerData.get(), _data.get() + offset, size);
                offset += size;

                if (pageManager != nullptr)
                    pageManager.reset();

                pageManager = std::make_unique<memory::PageManager<O>>();
                pageManager->deserialize(std::move(pageManagerData));

            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = sizeof(size_t) * 3; //page size, prunning and leafNodeAccess
            ans += sizeof(size_t) + indexName.size();
            ans += sizeof(size_t) + indexFolder.size();
            ans += sizeof(size_t) + buildFile.size();
            ans += sizeof(size_t) + experimentFile.size();
            ans += sizeof(size_t) + indexTypeMap[indexType].size();

            if (dataset != nullptr)
                ans += sizeof(size_t) + dataset->getSerializedSize();
            else
                ans += sizeof(size_t);

            if (distanceFunction != nullptr) {
                ans += sizeof(size_t) * 2 + distance::distanceTypeMap[distanceFunction->getDistanceType()].size() + distanceFunction->getSerializedSize();
            }
            else
                ans += sizeof(size_t);

            if (pivots != nullptr) {
                ans += sizeof(size_t) * 2 + pivots::PIVOT_TYPE2STR[pivots->getPivotType()].size() + pivots->getSerializedSize();
            }
            else
                ans += sizeof(size_t);

            if (pageManager != nullptr)
                ans += sizeof(size_t) + pageManager->getSerializedSize();
            else
                ans += sizeof(size_t);

            return ans;

        }

    };

    template <typename O, typename T>
    std::ostream& operator<<(std::ostream& os, const Index<O, T>& printable) {
        printable.print(os);
        return os;
    }
    
}

#endif //GERVLIB_INDEX_H
