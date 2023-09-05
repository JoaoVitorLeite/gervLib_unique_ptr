//
// Created by joaoleite on 8/30/23.
//

#ifndef GERVLIB_OMNIKDTREE_H
#define GERVLIB_OMNIKDTREE_H

#include "KdTree.h"

namespace gervLib::index::omni
{

    template <typename O, typename T>
    class OmniKdTree : public Index<O, T>
    {

    private:
        std::unique_ptr<dataset::Dataset<O, T>> mappedDataset;
        std::unique_ptr<kdtree::Node<O, T>> root;
        size_t numPerLeaf{}, numPivots{};
        bool storeLeafNode{}, storeDirectoryNode{}, useLAESA{};
        std::string serializedTree{};

    protected:
        std::string headerBuildFile() override
        {
            return "time,sys_time,user_time,distCount,iowrite,ioread";
        }

        std::string headerExperimentFile() override
        {
            return "expt_id,k,r,time,sys_time,user_time,distCount,prunning,iowrite,ioread";
        }

        double medianByDimension(std::unique_ptr<dataset::Dataset<O, T>>& dataset, size_t dim)
        {
            std::vector<double> vec(dataset->getCardinality());

            for (size_t i = 0; i < dataset->getCardinality(); i++)
                vec[i] = dataset->getElement(i).operator[](dim);

            std::sort(vec.begin(), vec.end());

            if(vec.size() % 2 == 0)
            {

                return (vec[vec.size()/2 - 1] + vec[vec.size()/2])/2.0;

            }
            else
            {

                return vec[vec.size()/2];

            }

        }

        std::pair<std::unique_ptr<std::vector<std::pair<double, double>>>, std::unique_ptr<std::vector<std::pair<double, double>>>> splitBounds(std::unique_ptr<std::vector<std::pair<double, double>>>& bound,
                                                                                                                                                double median, size_t pos)
        {
            std::unique_ptr<std::vector<std::pair<double, double>>> cp1 = std::make_unique<std::vector<std::pair<double, double>>>(*bound);
            std::unique_ptr<std::vector<std::pair<double, double>>> cp2 = std::make_unique<std::vector<std::pair<double, double>>>(*bound);
            std::pair<double, double> aux_cp1 = std::make_pair(bound->at(pos).first, std::nextafter(median, -std::numeric_limits<double>::infinity())),
                    aux_cp2 = std::make_pair(median, bound->at(pos).second);
            cp1->at(pos) = aux_cp1;
            cp2->at(pos) = aux_cp2;
            return std::make_pair(std::move(cp1), std::move(cp2));

        }
        void deleteRecursive(std::unique_ptr<kdtree::Node<O, T>> node)
        {
            if (node == nullptr)
                return;

            deleteRecursive(std::move(node->getLeft()));
            deleteRecursive(std::move(node->getRight()));
            node->clear();
            node.reset();

        }

        void clearRecursive(std::unique_ptr<kdtree::Node<O, T>>& node)
        {
            if (node == nullptr)
                return;

            clearRecursive(node->getLeft());
            clearRecursive(node->getRight());
            node->clear();

        }

        bool isEqualHelper(std::unique_ptr<kdtree::Node<O, T>>& node1, std::unique_ptr<kdtree::Node<O, T>>& node2)
        {
            if (node1 == nullptr && node2 == nullptr)
                return true;

            if ((node1 != nullptr && node2 == nullptr) || (node1 == nullptr && node2 != nullptr))
                return false;

            if (!node1->isEqual(node2))
                return false;

            return isEqualHelper(node1->getLeft(), node2->getLeft()) && isEqualHelper(node1->getRight(), node2->getRight());
        }

        std::string serializeTreeRecursive(std::unique_ptr<kdtree::Node<O, T>>& node)
        {
            if (node == nullptr)
                return "null ";

            if (node->getMemoryStatus() != index::MEMORY_STATUS::IN_DISK)
            {
                std::unique_ptr<u_char[]> data = node->serialize();
                this->pageManager->save(node->getNodeID(), std::move(data), node->getSerializedSize());
                node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);
            }

            std::string result = (node->isLeafNode() ? "L" : "D") + std::to_string(node->getNodeID()) + " 2 ";
            result += serializeTreeRecursive(node->getLeft());
            result += serializeTreeRecursive(node->getRight());

            return result;

        }

        std::unique_ptr<naryTree::NodeNAry> deserializeTreeRecursive(std::stringstream& ss)
        {
            std::string valStr;
            ss >> valStr;

            if (valStr == "null")
                return nullptr;

            size_t numChildren;
            ss >> numChildren;

            std::unique_ptr<naryTree::NodeNAry> node = std::make_unique<naryTree::NodeNAry>(valStr);

            for (size_t i = 0; i < numChildren; i++)
                node->children.push_back(deserializeTreeRecursive(ss));

            return node;
        }

        void buildTree(std::unique_ptr<kdtree::Node<O, T>>& node, std::unique_ptr<naryTree::NodeNAry>& aux)
        {

            if (aux == nullptr)
                return;

            if (aux->value[0] == 'L')
            {
                node = std::make_unique<kdtree::LeafNode<O, T>>();
                node->setNodeID(std::stoull(aux->value.substr(1)));

                if (!storeLeafNode) {
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_MEMORY);
                    std::unique_ptr<u_char[]> data = this->pageManager->load(node->getNodeID());
                    node->deserialize(std::move(data));
                }
                else
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);

                node->setLeft(nullptr);
                node->setRight(nullptr);
            }
            else
            {
                node = std::make_unique<kdtree::DirectoryNode<O, T>>();
                node->setNodeID(std::stoull(aux->value.substr(1)));

                if (!storeDirectoryNode) {
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_MEMORY);
                    std::unique_ptr<u_char[]> data = this->pageManager->load(node->getNodeID());
                    node->deserialize(std::move(data));
                }
                else
                    node->setMemoryStatus(index::MEMORY_STATUS::IN_DISK);

                node->setLeft(nullptr);
                node->setRight(nullptr);
            }

            buildTree(node->getLeft(), aux->children[0]);
            buildTree(node->getRight(), aux->children[1]);

        }

        bool isInterval(double infBound, double supBound, double test)
        {

            return ((test >= infBound) && (test <= supBound));

        }

        double minDist(dataset::BasicArrayObject<O, T>& query, std::unique_ptr<std::vector<std::pair<double, double>>>& bounds)
        {
            double limInfCase3 = -1.0;
            double limInfCase2 = -1.0;
            double answer = -1.0;
            bool within = true;

            for(size_t x = 0; x < bounds->size(); x++)
            {

                if(!isInterval(std::abs(bounds->at(x).first), std::abs(bounds->at(x).second), query.operator[](x)))
                {

                    within = false;

                    limInfCase3 = std::max(limInfCase3,

                                           std::min(

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).first)),

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).second))

                                           )

                    );
                }
                else
                {

                    limInfCase2 = std::min(limInfCase2,

                                           std::min(

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).first)),

                                                   std::abs(query.operator[](x) - std::abs(bounds->at(x).second))

                                           )

                    );

                }

            }

            if(within)
            {

                answer = 0.0;

            }
            else
            {

                if(limInfCase2 != -1.0)
                {

                    answer = limInfCase2;

                }
                else
                {

                    answer = limInfCase3;

                }

            }

            return answer;
        }

        double maxDist(dataset::BasicArrayObject<O, T>& query, std::unique_ptr<std::vector<std::pair<double, double>>>& bounds)
        {
            double answer = -1.0;

            answer = std::numeric_limits<double>::max();

            for(size_t x = 0; x < bounds->size(); x++)
            {

                if (std::numeric_limits<double>::max() -  std::abs(bounds->at(x).first) >= query.operator[](x))
                {

                    answer = std::min(answer, query.operator[](x) + std::abs(bounds->at(x).first));

                }

                if (std::numeric_limits<double>::max() -  std::abs(bounds->at(x).second) >= query.operator[](x))
                {

                    answer = std::min(answer, query.operator[](x) + std::abs(bounds->at(x).second));

                }

            }

            return answer;
        }

    public:
        OmniKdTree()
        {
            this->dataset = nullptr;
            this->mappedDataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->root = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::OMNIKDTREE;
            this->indexName = "OMNIKDTREE";
            this->indexFolder = "";
        }

        OmniKdTree(std::unique_ptr<dataset::Dataset<O, T>> _dataset,
                   std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> _df,
                   std::unique_ptr<pivots::Pivot<O, T>> _pivots, size_t _numPivots, size_t _numPerLeaf, size_t _pageSize = 0,
                   bool _storeDirectoryNode = false, bool _storeLeafNode = true, bool _useLAESA = false, std::string folder="")
        {

            this->dataset = std::move(_dataset);
            this->mappedDataset = nullptr;
            this->distanceFunction = std::move(_df);
            this->pivots = std::move(_pivots);
            this->root = nullptr;
            this->pageSize = _pageSize;
            this->numPivots = _numPivots;
            this->numPerLeaf = _numPerLeaf;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->storeLeafNode = _storeLeafNode;
            this->storeDirectoryNode = _storeDirectoryNode;
            this->useLAESA = _useLAESA;
            this->indexType = INDEX_TYPE::OMNIKDTREE;
            this->indexName = "OMNIKDTREE";

            if (!folder.empty())
                this->indexFolder = folder;

            this->generateIndexFiles(true, true);

            this->pageManager = std::make_unique<memory::PageManager<O>>("omnikd_page", this->indexFolder, this->pageSize);

            this->buildIndex();

        }

        explicit OmniKdTree(std::string _folder, std::string serializedFile = "")
        {
            this->dataset = nullptr;
            this->mappedDataset = nullptr;
            this->distanceFunction = nullptr;
            this->pivots = nullptr;
            this->pageManager = nullptr;
            this->root = nullptr;
            this->pageSize = 0;
            this->prunning = 0;
            this->leafNodeAccess = 0;
            this->numPerLeaf = 0;
            this->numPivots = 0;
            this->storeLeafNode = false;
            this->storeDirectoryNode = false;
            this->useLAESA = false;
            this->indexType = INDEX_TYPE::OMNIKDTREE;
            this->indexName = "OMNIKDTREE";
            this->indexFolder = _folder.empty() ? utils::generatePathByPrefix(configure::baseOutputPath, this->indexName) : _folder;

            if (serializedFile.empty())
                this->loadIndex();
            else
                this->loadIndex(serializedFile);

        }

        ~OmniKdTree() override
        {
            deleteRecursive(std::move(root));
        }

        std::unique_ptr<kdtree::Node<O, T>>& getRoot()
        {
            return root;
        }

        bool isEqual(std::unique_ptr<Index<O, T>>& other) override
        {
            if(!gervLib::index::Index<O, T>::isEqual(other))
                return false;

            auto* _other = dynamic_cast<OmniKdTree<O, T>*>(other.get());

            return isEqualHelper(this->root, _other->root);

        }

        void print(std::ostream& os) const override {

            std::stack<std::pair<kdtree::Node<O, T> *, size_t>> nodeStack;
            nodeStack.push(std::make_pair(root.get(), 0));

            os << "\n\n**********************************************************************************************************************************************************************************\n\n";
            os << "KdTree" << std::endl;
            os << "Number of pivots: " << numPivots << std::endl;
            os << "Number of objects per leaf: " << numPerLeaf << std::endl;
            os << "Store leaf node: " << (storeLeafNode ? "true" : "false") << std::endl;
            os << "Store directory node: " << (storeDirectoryNode ? "true" : "false") << std::endl;
            os << "Use LAESA: " << (useLAESA ? "true" : "false") << std::endl;

            while (!nodeStack.empty())
            {

                auto currentNode = nodeStack.top();
                nodeStack.pop();

                os << "**********************************************************************************************************************************************************************************\n\n";
                os << *currentNode.first << std::endl;

                if (currentNode.first->getLeft() != nullptr)
                    nodeStack.push(std::make_pair(currentNode.first->getLeft().get(), currentNode.second + 1));

                if (currentNode.first->getRight() != nullptr)
                    nodeStack.push(std::make_pair(currentNode.first->getRight().get(), currentNode.second + 1));

            }

        }

        void buildIndex() override {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            size_t currentNodeID = 0;
            std::queue<std::tuple<kdtree::Node<O, T> *, std::unique_ptr<dataset::Dataset<O, T>>, size_t>> nodeQueue;
            std::tuple<kdtree::Node<O, T>*, std::unique_ptr<dataset::Dataset<O, T>>, size_t> currentTuple;
            kdtree::Node<O, T> *currentNode;
            std::unique_ptr<dataset::Dataset<O, T>> currentDataset;
            size_t dPartition;
            double median;
            this->pivots->operator()(this->dataset, this->distanceFunction, this->numPivots);

            if (this->dataset->getCardinality() <= numPerLeaf)
                root = std::make_unique<kdtree::LeafNode<O, T>>();
            else
                root = std::make_unique<kdtree::DirectoryNode<O, T>>();

            root->setBoundsSize(this->numPivots);
            root->setNodeID(currentNodeID++);

            if (mappedDataset != nullptr)
            {
                mappedDataset->clear();
                mappedDataset.reset();
            }

            mappedDataset = std::make_unique<dataset::Dataset<O, T>>();
            mappedDataset->setSeed(this->dataset->getSeed());
            mappedDataset->setPath(this->dataset->getPath());
            mappedDataset->setDimensionality(this->dataset->getDimensionality());

            for (size_t i = 0; i < this->dataset->getCardinality(); i++)
            {
                dataset::BasicArrayObject<O, T> element = this->dataset->getElement(i);

                for (size_t j = 0; j < this->numPivots; j++)
                    element.operator[](j) = this->distanceFunction->operator()(this->dataset->getElement(i), this->pivots->getPivot(j));

                mappedDataset->insert(element);

            }

            nodeQueue.push(std::make_tuple(root.get(), std::move(this->mappedDataset), 0));

            while (!nodeQueue.empty()) {

                currentTuple = std::move(nodeQueue.front());
                nodeQueue.pop();

                currentNode = std::get<0>(currentTuple);
                currentDataset = std::move(std::get<1>(currentTuple));
                dPartition = std::get<2>(currentTuple);

                if (currentNode->isLeafNode())
                {

                    auto *leafNode = dynamic_cast<kdtree::LeafNode<O, T>*>(currentNode);

                    std::unique_ptr<dataset::Dataset<O, T>> auxDataset = std::make_unique<dataset::Dataset<O, T>>();
                    auxDataset->setSeed(currentDataset->getSeed());
                    auxDataset->setPath(currentDataset->getPath());
                    auxDataset->setDimensionality(currentDataset->getDimensionality());

                    for (size_t i = 0; i < currentDataset->getCardinality(); i++)
                    {
                        auxDataset->insert(this->dataset->getElement(currentDataset->getElement(i).getOID()));
                    }

                    if (useLAESA) {
                        std::filesystem::path leafIndexPath(this->indexFolder);
                        leafIndexPath /= "laesa_leafnode_" + std::to_string(currentNode->getNodeID());
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                this->distanceFunction->getDistanceType());
                        size_t oldDistCount = df->getDistanceCount();
                        std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                std::move(auxDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                this->numPivots,
                                leafIndexPath);
                        idx->getDistanceFunction()->setDistanceCount(idx->getDistanceFunction()->getDistanceCount() + oldDistCount);
                        leafNode->setIndex(std::move(idx));
                        currentDataset->clear();
                        currentDataset.reset();
                    }
                    else
                    {
                        leafNode->setDataset(std::move(auxDataset));
                        currentDataset->clear();
                        currentDataset.reset();
                    }
                }
                else
                {
                    median = medianByDimension(currentDataset, dPartition);

                    std::unique_ptr<dataset::Dataset<O, T>> leftDataset = std::make_unique<dataset::Dataset<O, T>>(), rightDataset = std::make_unique<dataset::Dataset<O, T>>();
                    leftDataset->setSeed(currentDataset->getSeed());
                    rightDataset->setSeed(currentDataset->getSeed());
                    leftDataset->setPath(currentDataset->getPath());
                    rightDataset->setPath(currentDataset->getPath());
                    leftDataset->setDimensionality(currentDataset->getDimensionality());
                    rightDataset->setDimensionality(currentDataset->getDimensionality());

                    for (size_t i = 0; i < currentDataset->getCardinality(); i++) {

                        if (currentDataset->getElement(i).operator[](dPartition) <= median)
                            leftDataset->insert(currentDataset->getElement(i));
                        else
                            rightDataset->insert(currentDataset->getElement(i));

                    }

                    std::pair<std::unique_ptr<std::vector<std::pair<double, double>>>, std::unique_ptr<std::vector<std::pair<double, double>>>> splitBoundaries = splitBounds(currentNode->getBoundary(), median, dPartition);

                    if (leftDataset->getCardinality() <= numPerLeaf)
                    {

                        std::unique_ptr<kdtree::LeafNode<O, T>> leafNode = std::make_unique<kdtree::LeafNode<O, T>>();
                        leafNode->setBoundary(std::move(splitBoundaries.first));
                        leafNode->setNodeID(currentNodeID++);

                        std::unique_ptr<dataset::Dataset<O, T>> auxDataset = std::make_unique<dataset::Dataset<O, T>>();
                        auxDataset->setSeed(currentDataset->getSeed());
                        auxDataset->setPath(currentDataset->getPath());
                        auxDataset->setDimensionality(currentDataset->getDimensionality());

                        for (size_t i = 0; i < leftDataset->getCardinality(); i++)
                        {
                            auxDataset->insert(this->dataset->getElement(leftDataset->getElement(i).getOID()));
                        }

                        if (useLAESA) {
                            std::filesystem::path leafIndexPath(this->indexFolder);
                            leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                    this->distanceFunction->getDistanceType());
                            size_t oldDistCount = df->getDistanceCount();
                            std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                    std::move(auxDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                    this->numPivots,
                                    leafIndexPath);
                            idx->getDistanceFunction()->setDistanceCount(idx->getDistanceFunction()->getDistanceCount() + oldDistCount);
                            leafNode->setIndex(std::move(idx));
                            leftDataset->clear();
                            leftDataset.reset();
                        }
                        else
                        {
                            leafNode->setDataset(std::move(auxDataset));
                            leftDataset->clear();
                            leftDataset.reset();
                        }

                        if (storeLeafNode)
                        {
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> leafData = leafNode->serialize();
                            this->pageManager->save(leafNode->getNodeID(), std::move(leafData), leafNode->getSerializedSize());
                            leafNode->clear();
                        }
                        else
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setLeft(std::move(leafNode));

                    }
                    else
                    {

                        std::unique_ptr<kdtree::DirectoryNode<O, T>> directoryNode = std::make_unique<kdtree::DirectoryNode<O, T>>();
                        directoryNode->setBoundary(std::move(splitBoundaries.first));
                        directoryNode->setNodeID(currentNodeID++);

                        if (this->storeDirectoryNode) {
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                            this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                            directoryNode->clear();
                        }
                        else
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setLeft(std::move(directoryNode));
                        nodeQueue.push(std::make_tuple(currentNode->getLeft().get(), std::move(leftDataset), (dPartition + 1) % currentDataset->getDimensionality()));

                    }

                    if (rightDataset->getCardinality() <= numPerLeaf)
                    {

                        std::unique_ptr<kdtree::LeafNode<O, T>> leafNode = std::make_unique<kdtree::LeafNode<O, T>>();
                        leafNode->setBoundary(std::move(splitBoundaries.second));
                        leafNode->setNodeID(currentNodeID++);

                        std::unique_ptr<dataset::Dataset<O, T>> auxDataset = std::make_unique<dataset::Dataset<O, T>>();
                        auxDataset->setSeed(currentDataset->getSeed());
                        auxDataset->setPath(currentDataset->getPath());
                        auxDataset->setDimensionality(currentDataset->getDimensionality());

                        for (size_t i = 0; i < rightDataset->getCardinality(); i++)
                        {
                            auxDataset->insert(this->dataset->getElement(rightDataset->getElement(i).getOID()));
                        }

                        if (useLAESA) {
                            std::filesystem::path leafIndexPath(this->indexFolder);
                            leafIndexPath /= "laesa_leafnode_" + std::to_string(leafNode->getNodeID());
                            std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>> df = distance::DistanceFactory<dataset::BasicArrayObject<O, T>>::createDistanceFunction(
                                    this->distanceFunction->getDistanceType());
                            size_t oldDistCount = df->getDistanceCount();
                            std::unique_ptr<Index<O, T>> idx = std::make_unique<index::LAESA<O, T>>(
                                    std::move(auxDataset), std::move(df), pivots::PivotFactory<O, T>::clone(this->pivots),
                                    this->numPivots, leafIndexPath);
                            idx->getDistanceFunction()->setDistanceCount(idx->getDistanceFunction()->getDistanceCount() + oldDistCount);
                            leafNode->setIndex(std::move(idx));
                            rightDataset->clear();
                            rightDataset.reset();
                        }
                        else
                        {
                            leafNode->setDataset(std::move(auxDataset));
                            rightDataset->clear();
                            rightDataset.reset();
                        }

                        if (storeLeafNode)
                        {
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> leafData = leafNode->serialize();
                            this->pageManager->save(leafNode->getNodeID(), std::move(leafData), leafNode->getSerializedSize());
                            leafNode->clear();
                        }
                        else
                            leafNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setRight(std::move(leafNode));

                    }
                    else
                    {

                        std::unique_ptr<kdtree::DirectoryNode<O, T>> directoryNode = std::make_unique<kdtree::DirectoryNode<O, T>>();
                        directoryNode->setBoundary(std::move(splitBoundaries.second));
                        directoryNode->setNodeID(currentNodeID++);

                        if (this->storeDirectoryNode) {
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                            std::unique_ptr<u_char[]> directoryData = directoryNode->serialize();
                            this->pageManager->save(directoryNode->getNodeID(), std::move(directoryData), directoryNode->getSerializedSize());
                            directoryNode->clear();
                        }
                        else
                            directoryNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);

                        currentNode->setRight(std::move(directoryNode));
                        nodeQueue.push(std::make_tuple(currentNode->getRight().get(), std::move(rightDataset), (dPartition + 1) % currentDataset->getDimensionality()));

                    }

                }

                if (currentNode->isLeafNode())
                {
                    if (storeLeafNode)
                    {
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> leafData = currentNode->serialize();
                        this->pageManager->save(currentNode->getNodeID(), std::move(leafData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);
                }
                else
                {
                    if (this->storeDirectoryNode) {
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_DISK);
                        std::unique_ptr<u_char[]> directoryData = currentNode->serialize();
                        this->pageManager->save(currentNode->getNodeID(), std::move(directoryData), currentNode->getSerializedSize());
                        currentNode->clear();
                    }
                    else
                        currentNode->setMemoryStatus(gervLib::index::MEMORY_STATUS::IN_MEMORY);
                }

                currentDataset->clear();
                currentDataset.reset();

            }

            this->dataset->clear();
            this->dataset.reset();

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

        void clear() override
        {
            clearRecursive(root);
        }

        std::vector<gervLib::query::ResultEntry<O>> kNNIncremental(gervLib::dataset::BasicArrayObject<O, T>& query, size_t k, bool saveResults) override
        {

            utils::Timer timer{};
            timer.start();
            this->distanceFunction->resetStatistics();
            this->prunning = 0;
            this->leafNodeAccess = 0;
            size_t ioW = configure::IOWrite, ioR = configure::IORead;
            std::priority_queue<query::Partition<kdtree::Node<O, T>*>, std::vector<query::Partition<kdtree::Node<O, T>*>>, std::greater<query::Partition<kdtree::Node<O, T>*>>> nodeQueue;
            std::priority_queue<query::ResultEntry<O>, std::vector<query::ResultEntry<O>>, std::greater<query::ResultEntry<O>>> elementQueue;
            query::Result<O> result;
            result.setMaxSize(k);
            query::Partition<kdtree::Node<O, T>*> currentPartition;
            kdtree::Node<O, T>* currentNode;
            kdtree::LeafNode<O, T>* currentLeafNode;
            LAESA<O, T>* laesa;
            double dist;
            dataset::BasicArrayObject<O, T> auxQuery = query;

            for (size_t i = 0; i < this->numPivots; i++)
            {
                auxQuery.operator[](i) = this->distanceFunction->operator()(query, this->pivots->getPivot(i));
            }

            nodeQueue.push(query::Partition<kdtree::Node<O, T>*>(root.get(), 0.0, std::numeric_limits<double>::max()));

            while (result.size() < k && !(nodeQueue.empty() && elementQueue.empty())) {

                if (elementQueue.empty())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->isLeafNode())
                    {

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                            currentNode->deserialize(std::move(nodeData));
                        }

                        currentLeafNode = (kdtree::LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->clear();
                        }

                    }
                    else
                    {

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getLeft()->getNodeID());
                            currentNode->getLeft()->deserialize(std::move(nodeData));
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getRight()->getNodeID());
                            currentNode->getRight()->deserialize(std::move(nodeData));
                        }

                        nodeQueue.push(query::Partition<kdtree::Node<O, T>*>(currentNode->getLeft().get(), minDist(auxQuery, currentNode->getLeft()->getBoundary()), maxDist(auxQuery, currentNode->getLeft()->getBoundary())));
                        nodeQueue.push(query::Partition<kdtree::Node<O, T>*>(currentNode->getRight().get(), minDist(auxQuery, currentNode->getRight()->getBoundary()), maxDist(auxQuery, currentNode->getRight()->getBoundary())));

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getLeft()->clear();
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getRight()->clear();
                        }

                    }

                }
                else if (!nodeQueue.empty() && nodeQueue.top().getMin() < elementQueue.top().getDistance())
                {

                    currentPartition = nodeQueue.top();
                    nodeQueue.pop();
                    currentNode = currentPartition.getElement();

                    if (currentNode->isLeafNode())
                    {

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getNodeID());
                            currentNode->deserialize(std::move(nodeData));
                        }

                        currentLeafNode = (kdtree::LeafNode<O, T>*) currentNode;
                        this->leafNodeAccess++;

                        if (currentLeafNode->getIndex() != nullptr)
                        {
                            laesa = (LAESA<O, T>*) currentLeafNode->getIndex().get();
                            std::vector<query::ResultEntry<O>> leafQuery = laesa->prunningQuery(query, k);

                            for (auto& entry : leafQuery)
                                elementQueue.push(entry);

                            this->prunning += currentLeafNode->getIndex()->getPrunning();

                        }
                        else
                        {

                            for (size_t i = 0; i < currentLeafNode->getDataset()->getCardinality(); i++)
                            {

                                dist = this->distanceFunction->operator()(query, currentLeafNode->getDataset()->getElement(i));
                                elementQueue.push(query::ResultEntry<O>(currentLeafNode->getDataset()->getElement(i).getOID(), dist));

                            }

                        }

                        if (currentNode->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->clear();
                        }

                    }
                    else
                    {

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getLeft()->getNodeID());
                            currentNode->getLeft()->deserialize(std::move(nodeData));
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            std::unique_ptr<u_char[]> nodeData = this->pageManager->load(currentNode->getRight()->getNodeID());
                            currentNode->getRight()->deserialize(std::move(nodeData));
                        }

                        nodeQueue.push(query::Partition<kdtree::Node<O, T>*>(currentNode->getLeft().get(), minDist(auxQuery, currentNode->getLeft()->getBoundary()), maxDist(auxQuery, currentNode->getLeft()->getBoundary())));
                        nodeQueue.push(query::Partition<kdtree::Node<O, T>*>(currentNode->getRight().get(), minDist(auxQuery, currentNode->getRight()->getBoundary()), maxDist(auxQuery, currentNode->getRight()->getBoundary())));

                        if (currentNode->getLeft()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getLeft()->clear();
                        }

                        if (currentNode->getRight()->getMemoryStatus() == MEMORY_STATUS::IN_DISK)
                        {
                            currentNode->getRight()->clear();
                        }

                    }

                }
                else
                {

                    result.push(elementQueue.top());
                    elementQueue.pop();

                }

            }

            std::vector<query::ResultEntry<O>> ans = result.getResults();
            std::reverse(ans.begin(), ans.end());

            std::string expt_id = utils::generateExperimentID();
            timer.stop();

            if (saveResults)
            {
                this->saveResultToFile(ans, query, "kNNIncremental", expt_id, {expt_id, std::to_string(k), "-1",
                                                                               std::to_string(timer.getElapsedTime()),
                                                                               std::to_string(timer.getElapsedTimeSystem()),
                                                                               std::to_string(timer.getElapsedTimeUser()),
                                                                               std::to_string(this->distanceFunction->getDistanceCount()),
                                                                               std::to_string(this->prunning),
                                                                               std::to_string(configure::IOWrite - ioW),
                                                                               std::to_string(configure::IORead - ioR)});
            }

            return ans;

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, sz;

            sz = gervLib::index::Index<O, T>::getSerializedSize();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = gervLib::index::Index<O, T>::serialize();
            memcpy(data.get() + offset, indexData.get(), sz);
            offset += sz;
            indexData.reset();

            memcpy(data.get() + offset, &numPerLeaf, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &numPivots, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storeDirectoryNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (storeLeafNode ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = (useLAESA ? 1 : 0);
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            sz = sizeof(size_t) + serializedTree.size();
            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, serializedTree.c_str(), serializedTree.size());
            offset += serializedTree.size();

            return data;
        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {
            size_t offset = 0, sz;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::unique_ptr<u_char[]> indexData = std::make_unique<u_char[]>(sz);
            memcpy(indexData.get(), _data.get() + offset, sz);
            offset += sz;
            gervLib::index::Index<O, T>::deserialize(std::move(indexData));

            memcpy(&numPerLeaf, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&numPivots, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeDirectoryNode = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            storeLeafNode = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            useLAESA = (sz == 1);

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(&aux[0], _data.get() + offset, sz);
            offset += sz;

            std::stringstream ss(aux);
            auto root_helper = deserializeTreeRecursive(ss);

            if (root != nullptr)
                root.reset();

            buildTree(root, root_helper);

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            this->serializedTree = serializeTreeRecursive(root);

            ans += sizeof(size_t) + gervLib::index::Index<O, T>::getSerializedSize();
            ans += sizeof(size_t) * 5;

            ans += sizeof(size_t) + serializedTree.size();

            return ans;
        }

    };

}

#endif //GERVLIB_OMNIKDTREE_H
