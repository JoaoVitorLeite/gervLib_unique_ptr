//
// Created by joaovictor on 18/08/23.
//

#ifndef GERVLIB_QUERY_H
#define GERVLIB_QUERY_H

#include <vector>
#include <queue>
#include <limits>
#include <functional>
#include <cstring>
#include "Serialize.h"
#include <format>
#include <iostream>

namespace gervLib::query
{

    template <typename T>
    concept HasSerializeSize = requires(T t) {
        {t.getSerializedSize()} -> std::convertible_to<size_t>;
    };

    template <typename T>
    concept IsSize_t = std::is_same_v<T, size_t>;

    template <typename O>
    requires HasSerializeSize<O> || IsSize_t<O>
    class ResultEntry: public gervLib::serialize::Serialize
    {

    private:
        O element;
        double distance{};

    public:
        ResultEntry() = default;

        ResultEntry(O element, double distance) : element(element), distance(distance) {}

        virtual ~ResultEntry() = default;

        O getElement() const { return element; }

        [[nodiscard]] double getDistance() const { return distance; }

        void setElement(O _element) { this->element = _element; }

        void setDistance(double _distance) { this->distance = _distance; }

        bool isEqual(const ResultEntry& e) const { return distance == e.distance; }

        bool operator<(const ResultEntry& e) const { return distance < e.distance; }

        bool operator>(const ResultEntry& e) const { return distance > e.distance; }

        bool operator==(const ResultEntry& e) const { return isEqual(e); }

        bool operator!=(const ResultEntry& e) const { return !isEqual(e); }

        bool operator<=(const ResultEntry& e) const { return distance <= e.distance; }

        bool operator>=(const ResultEntry& e) const { return distance >= e.distance; }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0;

            memcpy(data.get() + offset, &distance, sizeof(double));
            offset += sizeof(double);

            if constexpr (std::is_same_v<O, size_t>)
            {
                memcpy(data.get() + offset, &element, sizeof(O));
                offset += sizeof(O);
            }
            else
            {
                size_t size = element.getSerializedSize();
                memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> elementData = element.serialize();
                memcpy(data.get() + offset, elementData.get(), size);
                offset += size;
                elementData.reset();
            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0;

            memcpy(&distance, _data.get() + offset, sizeof(double));
            offset += sizeof(double);

            if constexpr (std::is_same_v<O, size_t>)
            {
                memcpy(&element, _data.get() + offset, sizeof(O));
                offset += sizeof(O);
            }
            else
            {
                size_t size;
                memcpy(&size, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> elementData(new u_char[size]);
                memcpy(elementData.get(), _data.get() + offset, size);
                offset += size;

                element.deserialize(std::move(elementData));
            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            if constexpr (std::is_same_v<O, size_t>)
                return sizeof(O) + sizeof(double);
            else
                return sizeof(size_t) + element.getSerializedSize() + sizeof(double);
        }

        friend std::ostream& operator<<(std::ostream& os, const ResultEntry& entry)
        {
            os << entry.distance << " <-> " << entry.element;
            return os;
        }

    };

    template <typename T>
    class CustomPriorityQueue {
    private:
        std::priority_queue<T, std::vector<T>, std::function<bool(const T&, const T&)>> pq;

    public:
        explicit CustomPriorityQueue(std::function<bool(const T&, const T&)> compareFunc) : pq(compareFunc) {}

        ~CustomPriorityQueue() {
            while (!pq.empty()) {
                pq.pop();
            }
        }

        void push(const T& value) {
            pq.push(value);
        }

        T top() {
            return pq.top();
        }

        void pop() {
            pq.pop();
        }

        bool empty() {
            return pq.empty();
        }

        size_t size() {
            return pq.size();
        }

        void clear() {
            while (!pq.empty()) {
                pq.pop();
            }
        }

        std::vector<T> dequeueInOrder()
        {

            std::vector<T> ans;
            std::priority_queue<T, std::vector<T>, std::function<bool(const T&, const T&)>> pqClone = pq;

            while(!pqClone.empty())
            {

                ans.push_back(pqClone.top());
                pqClone.pop();

            }

            return ans;

        }

        friend std::ostream& operator<<(std::ostream& os, const CustomPriorityQueue& _pq)
        {
            os << "Priority Queue:\n";
            std::priority_queue<T, std::vector<T>, std::function<bool(const T&, const T&)>> pqClone = _pq.pq;
            while(!pqClone.empty())
            {
                os << pqClone.top() << "\n";
                pqClone.pop();
            }
            return os;
        }

    };

    template <typename O>
    class Result : public gervLib::serialize::Serialize {

    private:
        size_t max_size;
        std::unique_ptr<CustomPriorityQueue<ResultEntry<O>>> pq;
        bool asc;

    public:
        explicit Result(bool _asc = false, size_t _max_size = std::numeric_limits<size_t>::max()) {

            asc = _asc;
            max_size = _max_size;

            if (asc)
                pq = std::make_unique<CustomPriorityQueue<ResultEntry<O>>>(
                        [](const ResultEntry<O> &e1, const ResultEntry<O> &e2) { return e1 > e2; });
            else
                pq = std::make_unique<CustomPriorityQueue<ResultEntry<O>>>(
                        [](const ResultEntry<O> &e1, const ResultEntry<O> &e2) { return e1 < e2; });

        }

        virtual ~Result() {
            while (!pq->empty()) {
                pq->pop();
            }

            pq.reset();
        }

        void push(ResultEntry<O> entry, double distance) {

            if (pq->empty() || pq->top().getDistance() > distance) {

                pq->push(entry);

                if (pq->size() > max_size)
                    pq->pop();

            }

        }

        void push(ResultEntry<O> entry) {

            pq->push(entry);

            if (pq->size() > max_size)
                pq->pop();

        }

        ResultEntry<O> top() { return pq->top(); }

        void pop() { pq->pop(); }

        std::vector<ResultEntry<O>> getResults() const { return pq->dequeueInOrder(); }

        [[nodiscard]] size_t size() const { return pq->size(); }

        [[nodiscard]] size_t getMaxSize() const { return max_size; }

        void setMaxSize(size_t _max_size) { this->max_size = _max_size; }

        void clear() {

            while (!pq->empty())
                pq->pop();

        }

        bool empty() { return pq->empty(); }

        bool isEqual(const Result<O> &other) {

            if (max_size != other.max_size)
                return false;

            if (asc != other.asc)
                return false;

            std::vector<ResultEntry<O>> results = getResults();
            std::vector<ResultEntry<O>> otherResults = other.getResults();

            if (results.size() != otherResults.size())
                return false;

            for (size_t i = 0; i < results.size(); i++) {

                if (results[i] != otherResults[i])
                    return false;

            }

            return true;

        }

        bool operator==(const Result<O> &other) { return isEqual(other); }

        bool operator!=(const Result<O> &other) { return !isEqual(other); }

        friend std::ostream &operator<<(std::ostream &os, const Result &result) {

            os << "Max Size: " << result.max_size << std::endl;
            os << "Size: " << result.size() << std::endl;
            os << *result.pq << std::endl;

            return os;

        }

        std::unique_ptr<u_char[]> serialize() override
        {

            std::vector<ResultEntry<O>> results = getResults();
            size_t offset = 0;
            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t size;

            memcpy(data.get() + offset, &max_size, sizeof(size_t));
            offset += sizeof(size_t);

            size = asc ? 1 : 0;
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            size = results.size();
            memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            for(auto& result : results)
            {
                std::unique_ptr<u_char[]> resultData = result.serialize();
                size = result.getSerializedSize();

                memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, resultData.get(), size);
                offset += size;
                resultData.reset();
            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0;
            size_t size, size_pq;

            memcpy(&max_size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);
            asc = size == 1;

            if(pq != nullptr)
            {
                pq->clear();
                pq.reset();
            }

            if (asc)
                pq = std::make_unique<CustomPriorityQueue<ResultEntry<O>>>(
                        [](const ResultEntry<O> &e1, const ResultEntry<O> &e2) { return e1 > e2; });
            else
                pq = std::make_unique<CustomPriorityQueue<ResultEntry<O>>>(
                        [](const ResultEntry<O> &e1, const ResultEntry<O> &e2) { return e1 < e2; });


            memcpy(&size_pq, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            for(size_t i = 0; i < size_pq; i++)
            {
                memcpy(&size, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> resultData(new u_char[size]);
                memcpy(resultData.get(), _data.get() + offset, size);
                offset += size;

                ResultEntry<O> result;
                result.deserialize(std::move(resultData));
                push(result);
            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            std::vector<ResultEntry<O>> results = getResults();

            size_t size = sizeof(size_t) * 3;

            for(auto& result : results)
                size += sizeof(size_t) + result.getSerializedSize();

            return size;

        }

    };


}

#endif //GERVLIB_QUERY_H
