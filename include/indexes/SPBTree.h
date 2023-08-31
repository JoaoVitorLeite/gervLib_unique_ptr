//
// Created by joaoleite on 8/31/23.
//

#ifndef GERVLIB_SPBTREE_H
#define GERVLIB_SPBTREE_H

#include "Index.h"
#include "IndexFactory.h"
#include "NAryTree.h"
#include "btree_multimap.hpp"
#include "HilbertCurve.h"
#include "EquiDepth.h"

namespace gervLib::index::spb
{

    template <typename Key, typename Value>
    struct btree_spb_traits {
        //! If true, the tree will self verify its invariants after each insert() or
        //! erase(). The header must have been compiled with TLX_BTREE_DEBUG
        //! defined.
        static const bool self_verify = false;

        //! If true, the tree will print out debug information and a tree dump
        //! during insert() or erase() operation. The header must have been
        //! compiled with TLX_BTREE_DEBUG defined and key_type must be std::ostream
        //! printable.
#ifdef TLX_BTREE_DEBUG
        static const bool debug = true;
#else
        static const bool debug = false;
#endif

        //! Number of slots in each leaf of the tree. Estimated so that each node
        //! has a size of about 256 bytes.
//        static const int leaf_slots =
//                TLX_BTREE_MAX(8, 256 / (sizeof(Value)));
        static const int leaf_slots = 5;

        //! Number of slots in each inner node of the tree. Estimated so that each
        //! node has a size of about 256 bytes.
        static const int inner_slots =
                TLX_BTREE_MAX(8, 256 / (sizeof(Key) + sizeof(void*)));

        //! As of stx-btree-0.9, the code does linear search in find_lower() and
        //! find_upper() instead of binary_search, unless the node size is larger
        //! than this threshold. See notes at
        //! http://panthema.net/2013/0504-STX-B+Tree-Binary-vs-Linear-Search
        static const size_t binsearch_threshold = 256;
    };


    template <typename O, typename T, typename D>
    class SPBTree : public Index<O, T>
    {

    private:
        typedef tlx::btree_multimap<D, dataset::BasicArrayObject<O, T>, O, T, std::less<D>, btree_spb_traits<D, size_t>> btree_type;
        std::unique_ptr<equidepth::EquiDepth<T>> equiDepth;
        std::unique_ptr<hilbert::HilbertCurve<D>> hilbertCurve;
        size_t numPerLeaf{}, numPivots{}, num_bins{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::unique_ptr<btree_type> bptree;

    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount,iowrite,ioread";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning,iowrite,ioread";
        }

    public:
        SPBTree()
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->indexFolder = "";
            this->bptree = nullptr;
            this->equiDepth = nullptr;
            this->hilbertCurve = nullptr;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->num_bins = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
        }

        SPBTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
                std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
                std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _num_bins, size_t _numPerLeaf, size_t _pageSize = 0,
                bool _storeDirectoryNode = false, bool _storeLeafNode = false, bool _useLAESA = true, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->numPivots = _numPivots;
            this->numPerLeaf = _numPerLeaf;
            this->pageSize = _pageSize;
            this->num_bins = _num_bins;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->storeLeafNode = _storeLeafNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";
            this->bptree = nullptr;

            if (numPivots > 7)
                if constexpr (!configure::is_mpz_class_v<D>)
                    throw std::runtime_error("SPBTree: Number of pivots must be less than 8");

            unsigned long long p = log2(num_bins-1)+1;

            if constexpr (!configure::is_mpz_class_v<D>)
            {
                this->hilbertCurve = std::make_unique<hilbert::HilbertCurve_ull>(p, this->numPivots);
            }
            else
            {
                this->hilbertCurve = std::make_unique<hilbert::HilbertCurve_mpz>(mpz_class(std::to_string(p)), mpz_class(std::to_string(this->numPivots)));
            }

            this->equiDepth = std::make_unique<equidepth::EquiDepth<T>>(this->num_bins, this->numPivots);

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("spb_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit SPBTree(std::string _folder, std::string serializedFile = "")
        {
            this->dataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->bptree = nullptr;
            this->equiDepth = nullptr;
            this->hilbertCurve = nullptr;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->num_bins = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::SPBTREE;
            this->indexName = "SPBTREE";

            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~SPBTree() override = default;

        //TODO implement delete, clear, print, isEqual, buildIndex, kNN, kNNIncremental, serialize, deserialize, getSerializedSize

        void buildIndex() override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);
            this->bptree = std::make_unique<btree_type>();

            std::vector<std::vector<double>> pivotMapping = std::vector<std::vector<double>>(this->dataset->getCardinality(), std::vector<double>(this->numPivots, 0.0));
            std::vector<std::vector<D>> discretePivotMapping = std::vector<std::vector<D>>(this->dataset->getCardinality(), std::vector<D>(this->numPivots, 0));
            std::vector<D> keys(this->dataset->getCardinality());

            std::ofstream file_pivot_mapping(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_pivot_mapping.txt");
            std::ofstream file_disc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_disc.txt");
            std::ofstream file_sfc(this->indexFolder + std::filesystem::path::preferred_separator + "bulk_load_sfc.txt");

            for(size_t i = 0; i < this->dataset->getCardinality(); i++) {

                for (size_t j = 0; j < numPivots; j++) {

                    pivotMapping[i][j] = this->distanceFunction->getDistance(this->dataset->getElement(i), this->pivots->getPivot(j));
                    file_pivot_mapping << pivotMapping[i][j] << " ";

                }

                file_pivot_mapping << "\n";

            }

            this->equiDepth->build(pivotMapping);
            this->equiDepth->saveToFile(this->indexFolder);
            file_pivot_mapping.close();

            for(size_t i = 0; i < this->dataset->getCardinality(); i++)
            {

                for(size_t j = 0; j < numPivots; j++)
                {

                    discretePivotMapping[i][j] = this->equiDepth->getBin(j, pivotMapping[i][j]);
                    file_disc << discretePivotMapping[i][j] << " ";

                }

                file_disc << "\n";

                keys[i] = this->hilbertCurve->distance_from_point(discretePivotMapping[i]);
                file_sfc << keys[i] << "\n";

            }

            pivotMapping.clear();
            discretePivotMapping.clear();
            this->equiDepth->clear();
            file_disc.close();
            file_sfc.close();

            std::vector<std::pair<D, dataset::BasicArrayObject<O, T>>> insertValues(keys.size());

            for(size_t i = 0; i < keys.size(); i++)
            {

                insertValues[i] = std::make_pair(keys[i], this->dataset->getElement(i));

            }

            keys.clear();

            std::sort(insertValues.begin(), insertValues.end());
            bptree->bulk_load(insertValues.begin(), insertValues.end());
            insertValues.clear();
            bptree->setVariables(std::move(this->pivots), std::move(this->distanceFunction), std::move(this->hilbertCurve), std::move(this->pageManager),
                                 this->numPivots, this->indexFolder, this->storeDirectoryNode, this->storeLeafNode, this->useLAESA);
            bptree->build_MBR();
            bptree->test();











            timer.stop();

            std::ofstream buildFile(this->buildFile, std::ios::app);
            buildFile << timer.getElapsedTime()
                      << ","
                      << timer.getElapsedTimeSystem()
                      << ","
                      << timer.getElapsedTimeUser()
                      << ","
                      << this->distanceFunction->getDistanceCount()
                      << ","
                      << std::to_string(configure::IOWrite - ioW)
                      << ","
                      << std::to_string(configure::IORead - ioR) << std::endl;
            buildFile.close();

        }

        std::vector<gervLib::query::ResultEntry<O>> kNN(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("Not implemented yet");
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {
            throw std::runtime_error("Not implemented yet");
        }


    };



}

#endif //GERVLIB_SPBTREE_H
