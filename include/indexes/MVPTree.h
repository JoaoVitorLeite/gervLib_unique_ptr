//
// Created by joaovictor on 25/08/23.
//

#ifndef GERVLIB_MVPTREE_H
#define GERVLIB_MVPTREE_H

#include "Index.h"

namespace gervLib::index::mvptree
{

    template <typename O, typename T>
    class Node : public serialize::Serialize
    {

    protected:
        std::unique_ptr<std::vector<dataset::BasicArrayObject<O, T>>> pivots;
        std::vector<std::unique_ptr<Node<O, T>>> childrens;
        size_t nodeID{};
        std::unique_ptr<std::vector<std::vector<double>>> splits;
        MEMORY_STATUS memoryStatus = MEMORY_STATUS::NONE;

    public:
        Node() = default;

        explicit Node(size_t fo)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
        }

        Node(size_t lpn, size_t fo)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
            pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(lpn);
        }

        Node(size_t lpn, size_t fo, size_t ns)
        {
            childrens = std::vector<std::unique_ptr<Node<O, T>>>(fo);
            pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(lpn);
            splits = std::make_unique<std::vector<std::vector<double>>>(lpn, std::vector<double>(ns, -1.0));
        }

        virtual ~Node()
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }

            if (splits != nullptr)
            {
                splits->clear();
                splits.reset();
            }
        }

        void setChildren(size_t index, std::unique_ptr<Node<O, T>> node)
        {
            utils::check_range(0, childrens.size() - 1, index, "Node::setChildren(): Index out of range.");
            childrens[index] = std::move(node);
        }

        void setPivot(size_t index, dataset::BasicArrayObject<O, T> pivot)
        {

            if (pivots == nullptr)
                throw std::runtime_error("Node::setPivot(): Pivots vector is null.");

            utils::check_range(0, childrens.size() - 1, index, "Node::setChildren(): Index out of range.");
            pivots->at(index) = pivot;
        }

        void setSplit(size_t r, size_t c, double value)
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::setSplit(): Splits vector is null.");

            utils::check_range(0, splits->size() - 1, r, "Node::setSplit(): Index out of range.");
            utils::check_range(0, splits->at(r).size() - 1, c, "Node::setSplit(): Index out of range.");

            splits->at(r)[c] = value;

        }

        void setNodeID(size_t id) { nodeID = id; }

        std::unique_ptr<Node<O, T>>& getChildren(size_t index)
        {
            utils::check_range(0, childrens.size() - 1, index, "Node::getChildren(): Index out of range.");
            return childrens[index];
        }

        dataset::BasicArrayObject<O, T>& getPivot(size_t index)
        {
            if (pivots == nullptr)
                throw std::runtime_error("Node::getPivot(): Pivots vector is null.");

            utils::check_range(0, pivots->size() - 1, index, "Node::getPivot(): Index out of range.");
            return pivots->at(index);
        }

        double getSplit(size_t r, size_t c)
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::getSplit(): Splits vector is null.");

            utils::check_range(0, splits->size() - 1, r, "Node::getSplit(): Index out of range.");
            utils::check_range(0, splits->at(r).size() - 1, c, "Node::getSplit(): Index out of range.");

            return splits->at(r)[c];
        }

        size_t getNumberOfChildrens() { return childrens.size(); }

        size_t getNumberOfPivots()
        {
            if (pivots == nullptr)
                throw std::runtime_error("Node::getNumberOfPivots(): Pivots vector is null.");

            return pivots->size();
        }

        size_t getNumberOfSplits()
        {
            if (splits == nullptr)
                throw std::runtime_error("Node::getNumberOfSplits(): Splits vector is null.");

            return splits->at(0).size();
        }

        size_t getNodeID() { return nodeID; }

        bool operator==(Node<O, T> *other)
        {
            return isEqual(other);
        }

        bool operator!=(Node<O, T> *other)
        {
            return !isEqual(other);
        }

        //virtual methods

        virtual void clear()
        {
            if (pivots != nullptr)
            {
                pivots->clear();
                pivots.reset();
            }

            if (splits != nullptr)
            {
                splits->clear();
                splits.reset();
            }
        }

        virtual bool isLeaf()
        {
            return false;
        }

        virtual bool isEqual(std::unique_ptr<Node<O, T>>& other)
        {
            if (memoryStatus == MEMORY_STATUS::IN_DISK || other->memoryStatus == MEMORY_STATUS::IN_DISK)
                return nodeID == other->nodeID;
            else
            {

                if (((pivots == nullptr && other->pivots != nullptr) || (pivots != nullptr && other->pivots == nullptr)) || (pivots->size() != other->pivots->size()))
                    return false;
                else
                {
                    for (size_t i = 0; i < pivots->size(); i++)
                    {
                        if (pivots->at(i) != other->pivots->at(i))
                            return false;
                    }
                }

                if (((splits == nullptr && other->splits != nullptr) || (splits != nullptr && other->splits == nullptr)) || (splits->size() != other->splits->size() || splits->at(0).size() != other->splits->at(0).size()))
                    return false;
                else
                {
                    for (size_t i = 0; i < splits->size(); i++)
                    {
                        if (splits->at(i) != other->splits->at(i))
                            return false;
                    }
                }

                return true;

            }
        }

        virtual void print(std::ostream& os) const
        {
            os << "Node type: Node" << std::endl;
            os << "Node ID: " << nodeID << std::endl;

            if (pivots != nullptr)
            {
                os << "Pivots: " << std::endl;
                for (size_t i = 0; i < pivots->size(); i++)
                    os << pivots->at(i) << std::endl;
            }
            else
                os << "Pivots: null" << std::endl;

            if (splits != nullptr)
            {
                os << "Splits: " << std::endl;
                for (size_t i = 0; i < splits->size(); i++)
                {
                    for (size_t j = 0; j < splits->at(i).size(); j++)
                        os << splits->at(i)[j] << " ";
                    os << std::endl;
                }
            }
            else
                os << "Splits: null" << std::endl;

            os << "Memory Status: " << gervLib::index::memoryStatusMap[memoryStatus] << std::endl;

        }

        std::unique_ptr<u_char[]> serialize() override
        {
            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, sz;

            memcpy(data.get() + offset, &nodeID, sizeof(size_t));
            offset += sizeof(size_t);

            if (pivots != nullptr)
            {
                sz = pivots->size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                for (auto pvt : *pivots)
                {
                    sz = pvt.getSerializedSize();
                    memcpy(data.get() + offset, &sz, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> pvtData = pvt.serialize();
                    memcpy(data.get() + offset, pvtData.get(), sz);
                    offset += sz;
                    pvtData.reset();
                }
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            if (splits != nullptr)
            {
                sz = splits->size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                sz = splits->at(0).size();
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                for (auto split : *splits)
                {
                    memcpy(data.get() + offset, split.data(), sizeof(double) * split.size());
                    offset += sizeof(double) * split.size();
                }
            }
            else
            {
                sz = 0;
                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, &sz, sizeof(size_t));
                offset += sizeof(size_t);
            }

            std::string aux = gervLib::index::memoryStatusMap[memoryStatus];
            sz = aux.size();

            memcpy(data.get() + offset, &sz, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, aux.c_str(), sz);
            offset += sz;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, sz;

            memcpy(&nodeID, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            size_t sz_loop;
            memcpy(&sz_loop, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz_loop != 0)
            {

                if (pivots != nullptr)
                {
                    pivots->clear();
                    pivots.reset();
                }

                pivots = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(sz_loop);

                for (size_t i = 0; i < sz_loop; i++)
                {
                    memcpy(&sz, _data.get() + offset, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> pvtData(new u_char[sz]);
                    memcpy(pvtData.get(), _data.get() + offset, sz);
                    offset += sz;

                    pivots->at(i).deserialize(std::move(pvtData));
                }
            }
            else
                pivots = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            size_t sz2;
            memcpy(&sz2, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (sz != 0 && sz2 != 0)
            {

                if (splits != nullptr)
                {
                    splits->clear();
                    splits.reset();
                }

                splits = std::make_unique<std::vector<std::vector<double>>>(sz, std::vector<double>(sz2));

                for (size_t i = 0; i < sz; i++)
                {
                    memcpy(splits->at(i).data(), _data.get() + offset, sizeof(double) * sz2);
                    offset += sizeof(double) * sz2;
                }
            }
            else
                splits = nullptr;

            memcpy(&sz, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            std::string aux;
            aux.resize(sz);

            memcpy(aux.data(), _data.get() + offset, sz);
            offset += sz;

            memoryStatus = gervLib::index::memoryStatusMapReverse[aux];

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t ans = 0;

            ans += sizeof(size_t); //nodeID

            ans += sizeof(size_t); //pivots->size()
            if (pivots != nullptr)
            {

                for (auto pvt : *pivots)
                    ans += sizeof(size_t) + pvt.getSerializedSize();

            }

            ans += sizeof(size_t); //splits->size()
            ans += sizeof(size_t); //splits->at(0).size()
            if (splits != nullptr)
            {
                ans += sizeof(double) * splits->size() * splits->at(0).size();
            }

            ans += sizeof(size_t) + memoryStatusMap[memoryStatus].size(); //memoryStatus

            return ans;

        }

    };

}

#endif //GERVLIB_MVPTREE_H
