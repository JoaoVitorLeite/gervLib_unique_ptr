//
// Created by joaovictor on 15/08/23.
//

#ifndef GERVLIB_BASICARRAYOBJECT_H
#define GERVLIB_BASICARRAYOBJECT_H

#include "Serialize.h"
#include "Utils.h"
#include "Configure.h"
#include <vector>
//#include <format>
#include <numeric>
#include <cstring>

namespace gervLib::dataset
{
    template <typename O, typename T>
    class BasicArrayObject : public serialize::Serialize
    {

    private:
        std::unique_ptr<std::vector<T>> data;
        O oid;

    private:
        std::unique_ptr<u_char[]> serializeNumeric()
        {

            size_t size = getSerializedSize();
            std::unique_ptr<u_char[]> _data = std::make_unique<u_char[]>(size);
            size_t offset = 0;

            std::memcpy(_data.get() + offset, &oid, sizeof(O));
            offset += sizeof(O);

            size = this->size();
            std::memcpy(_data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            for(size_t x = 0; x < size; x++)
            {

                std::memcpy(_data.get() + offset, &data->at(x), sizeof(T));
                offset += sizeof(T);

            }

            return _data;

        }

        void deserializeNumeric(std::unique_ptr<u_char[]> _data)
        {

            size_t offset = 0;

            std::memcpy(&oid, _data.get() + offset, sizeof(O));
            offset += sizeof(O);

            size_t size;
            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
                return;

            this->clear();
            data = std::make_unique<std::vector<T>>(size);

            for(size_t x = 0; x < size; x++)
            {

                std::memcpy(&data->at(x), _data.get() + offset, sizeof(T));
                offset += sizeof(T);

            }

            _data.reset();

        }

        std::unique_ptr<u_char[]> serializeVecChar()
        {

            size_t size = getSerializedSize();
            std::unique_ptr<u_char[]> _data = std::make_unique<u_char[]>(size);
            size_t offset = 0;

            std::memcpy(_data.get() + offset, &oid, sizeof(O));
            offset += sizeof(O);

            size = this->size();
            std::memcpy(_data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            size_t auxSize;

            for(size_t x = 0; x < size; x++)
            {

                auxSize = data->at(x).size();
                std::memcpy(_data.get() + offset, &auxSize, sizeof(size_t));
                offset += sizeof(size_t);

                std::memcpy(_data.get() + offset, data->at(x).data(), auxSize * sizeof(char));
                offset += auxSize * sizeof(char);

            }

            return _data;

        }

        void deserializeVecChar(std::unique_ptr<u_char[]> _data)
        {

            size_t offset = 0;

            std::memcpy(&oid, _data.get() + offset, sizeof(O));
            offset += sizeof(O);

            size_t size;
            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
                return;

            this->clear();
            data = std::make_unique<std::vector<T>>(size);

            size_t auxSize;

            for(size_t x = 0; x < size; x++)
            {

                std::memcpy(&auxSize, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                data->at(x).resize(auxSize);
                std::memcpy(data->at(x).data(), _data.get() + offset, auxSize * sizeof(char));
                offset += auxSize * sizeof(char);

            }

            _data.reset();

        }

    public:
        //Constructors and destructor

        BasicArrayObject() : data(std::make_unique<std::vector<T>>()), oid(std::numeric_limits<O>::max()) {}

        BasicArrayObject(const O oid, const size_t size): data(std::make_unique<std::vector<T>>(size)), oid(oid) {}

        BasicArrayObject(const O oid, std::vector<T> data): data(std::make_unique<std::vector<T>>(data)), oid(oid) {}

        BasicArrayObject(const BasicArrayObject& other)
        {
            oid = other.oid;
            if (other.data != nullptr)
                data = std::make_unique<std::vector<T>>(*other.data);
        }

        virtual ~BasicArrayObject()
        {
            if (data != nullptr)
                data.reset();
        }

        // Virtual methods

        std::unique_ptr<u_char[]> serialize() override
        {

            if constexpr (std::is_same<T, double>::value || std::is_same<T, float>::value)
            {

                return serializeNumeric();

            }
            else if(std::is_same<T, std::vector<char>>::value)
            {

                return serializeVecChar();

            }
            else
            {

                throw std::runtime_error("BasicArrayObject::serialize: Type not supported");

            }

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            if constexpr (std::is_same<T, double>::value || std::is_same<T, float>::value)
            {

                deserializeNumeric(std::move(_data));

            }
            else if(std::is_same<T, std::vector<char>>::value)
            {

                deserializeVecChar(std::move(_data));

            }
            else
            {

                throw std::runtime_error("BasicArrayObject::deserialize: Type not supported");

            }

        }

        size_t getSerializedSize() override
        {
            if constexpr (std::is_same<T, double>::value || std::is_same<T, float>::value)
            {

                return sizeof(O) + sizeof(size_t) + (sizeof(T) * this->size());

            }
            else if(std::is_same<T, std::vector<char>>::value)
            {

                size_t total = sizeof(O) + sizeof(size_t);

                for(size_t x = 0; x < this->size(); x++)
                {

                    total += sizeof(size_t) + data->at(x).size() * sizeof(char);

                }

                return total;

            }
            else
            {

                throw std::runtime_error("BasicArrayObject::getSerializedSize: Type not supported");

            }
        }

        // Setters

        void setOID(const O _oid)
        {
            this->oid = _oid;
        }

        void setData(const std::vector<T> _data)
        {
            if (data == nullptr)
                data = std::make_unique<std::vector<T>>(_data);
            else {
                data->clear();
                data->insert(data->begin(), _data.begin(), _data.end());
            }
        }

        // Getters

        O getOID() const
        {
            return oid;
        }

        std::vector<T> getData() const
        {
            if (data == nullptr)
                return std::vector<T>();
            else
                return *data;
        }

        // Operators

        T operator[](const size_t index) const
        {
            if (data == nullptr)
                throw std::runtime_error("BasicArrayObject::operator[]: data is null");

            utils::check_range(0, size()-1, index, "BasicArrayObject::operator[]");
            return (*data)[index];
        }

        T& operator[](const size_t index)
        {
            if (data == nullptr)
                throw std::runtime_error("BasicArrayObject::operator[]: data is null");

            utils::check_range(0, size()-1, index, "BasicArrayObject::operator[]");
            return (*data)[index];
        }

        friend std::ostream& operator<<(std::ostream& os, const BasicArrayObject& obj)
        {

//            os << std::format("| {1:^{0}} || ", configure::oid_size, obj.oid);
//
//            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
//            {
//
//                for(size_t x = 0; x < obj.size(); x++)
//                {
//
//                    os << std::format("{1:^{0}} | ", configure::double_size, (*obj.data)[x]);
//
//                }
//
//            }
//            else
//            {
//
//                size_t maxSize = 0;
//
//                for(size_t x = 0; x < obj.size(); x++)
//                {
//
//                    if(obj.data->at(x).size() > maxSize)
//                    {
//
//                        maxSize = obj.data->at(x).size();
//
//                    }
//
//                }
//
//                maxSize += 2;
//                maxSize = std::max(maxSize, configure::str_size);
//                std::vector<char> aux;
//
//                for(size_t x = 0; x < obj.size(); x++)
//                {
//
//                    aux = (*obj.data)[x];
//
//                    os << std::format("{1:^{0}} | ", maxSize, std::accumulate(aux.begin(), aux.end(), std::string()));
//
//                    aux.clear();
//
//                }
//
//            }

            return os;

        }

        BasicArrayObject& operator=(const BasicArrayObject& other)
        {
            if (this == &other)
                return *this;

            oid = other.oid;
            if (other.data != nullptr)
                data = std::make_unique<std::vector<T>>(*other.data);
            else
                data.reset();

            return *this;
        }

        bool operator==(const BasicArrayObject& other) const
        {
            return isEqual(other);
        }

        bool operator!=(const BasicArrayObject& other) const
        {
            return !isEqual(other);
        }

        bool operator<(const BasicArrayObject& other) const
        {
            return oid < other.oid;
        }

        bool operator>(const BasicArrayObject& other) const
        {
            return oid > other.oid;
        }

        bool operator<=(const BasicArrayObject& other) const
        {
            return oid <= other.oid;
        }

        bool operator>=(const BasicArrayObject& other) const
        {
            return oid >= other.oid;
        }

        //Methods

        void clear()
        {
            if (data != nullptr)
                data.reset();
        }

        [[nodiscard]] size_t size() const
        {
            if (data == nullptr)
                return 0;
            else
                return data->size();
        }

        void set(size_t pos, T value)
        {
            if (data == nullptr)
                throw std::runtime_error("BasicArrayObject::set: data is null");

            gervLib::utils::check_range(0, size()-1, pos, "BasicArrayObject::set");
            data->at(pos) = value;
        }

        void set(T value)
        {
            if (data == nullptr)
                data = std::make_unique<std::vector<T>>(1, value);
            else
                data->push_back(value);
        }

        std::unique_ptr<BasicArrayObject<O, T>> clone() const
        {
            if (data == nullptr)
                return std::make_unique<BasicArrayObject<O, T>>(oid, 0);
            else
                return std::make_unique<BasicArrayObject<O, T>>(oid, *data);
        }

        bool isEqual(const BasicArrayObject& other) const
        {
            if (size() != other.size())
                return false;

            if (oid != other.oid)
                return false;

            if (data != nullptr)
                if (!std::equal(data->begin(), data->end(), other.data->begin()))
                    return false;

            return true;
        }

    };

}

#endif //GERVLIB_BASICARRAYOBJECT_H
