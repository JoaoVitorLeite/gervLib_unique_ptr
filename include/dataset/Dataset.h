//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_DATASET_H
#define GERVLIB_DATASET_H

#include "BasicArrayObject.h"
#include <regex>
#include <fstream>
#include <sstream>
#include <cmath>

namespace gervLib::dataset
{

    template <typename O, typename T>
    class Dataset : public serialize::Serialize
    {

    private:
        std::unique_ptr<std::vector<dataset::BasicArrayObject<O, T>>> dataset;
        size_t dimensionality, seed;
        std::string path;

    private:
        void read(size_t cardinality, size_t _dimensionality, const std::string& pattern = "")
        {

            if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
            {

                readNumeric(cardinality, _dimensionality, pattern);

            }
            else if (std::is_same_v<T, std::vector<char>>)
            {

                readText(cardinality, _dimensionality, pattern);

            }
            else
                throw std::runtime_error("Dataset::read: Type not supported");

        }

        template<typename W>
        W convertFromString(const std::string& str) {
            std::istringstream iss(str);
            W value;
            iss >> value;
            return value;
        }

        template<typename U>
        std::vector<U> split(const std::string& str, const std::string& delim)
        {

            std::regex rgx(delim);
            std::sregex_token_iterator iter(str.begin(),str.end(),rgx,-1);
            std::sregex_token_iterator end;
            std::vector<U> v;
            U d;

            for(; iter != end; ++iter)
            {

                d = convertFromString<U>(*iter);
                v.push_back(d);

            }

            return v;

        }

        std::vector<std::string> splitString(const std::string& str, const std::string& delim)
        {

            std::regex rgx(delim);
            std::sregex_token_iterator iter(str.begin(),str.end(),rgx,-1);
            std::sregex_token_iterator end;
            std::vector<std::string> v;

            for(; iter != end; ++iter)
            {

                v.push_back(*iter);

            }

            return v;

        }

        void readNumeric(size_t cardinality, size_t _dimensionality, const std::string& pattern = "")
        {

            if((cardinality != 0) && (_dimensionality != 0))
            {

                std::ifstream file(path);
                T readin;

                if (file.is_open())
                {

                    for (size_t x = 0; x < cardinality; x++)
                    {

                        dataset::BasicArrayObject<O, T> aux = dataset::BasicArrayObject<O, T>();
                        aux.setOID(x);
                        for (size_t y = 0; y < _dimensionality; y++)
                        {

                            file >> readin;
                            aux.set(readin);

                        }
                        dataset->push_back(aux);

                    }

                }
                else
                    throw std::runtime_error("Dataset::readNumeric: File not found");


                file.close();


            }
            else if(!pattern.empty())
            {

                std::ifstream file;
                std::string line;
                std::vector<T> aux;

                file.open(path);

                if (file.is_open()){

                    size_t oid = 1;

                    while(std::getline(file, line) && oid)
                    {

                        aux = split<T>(line, pattern);

                        dataset::BasicArrayObject<O, T> inst = dataset::BasicArrayObject<O, T>(oid++-1, aux);
                        dataset->push_back(inst);

                    }

                    setDimensionality(aux.size());

                }
                else
                    throw std::runtime_error("Dataset::readNumeric: File not found");

                file.close();

            }
            else
                throw std::runtime_error("Dataset::readNumeric: Invalid parameters");

        }

        void readText(size_t cardinality, size_t _dimensionality, const std::string& pattern = "")
        {

            if((cardinality != 0) && (_dimensionality != 0))
            {

                std::ifstream file;

                file.open(path);

                if (file.is_open())
                {

                    for(size_t x = 0; x < cardinality; x++)
                    {

                        dataset::BasicArrayObject<O, T> aux;
                        aux.setOID(x);

                        for(size_t y = 0; y < _dimensionality; y++)
                        {

                            std::string s;
                            file >> s;

                            T charVec(s.begin(), s.end());
                            aux.set(charVec);

                        }

                        dataset->push_back(aux);

                    }

                }
                else
                    throw std::runtime_error("Dataset::readText: File not found");


            }
            else if(!pattern.empty())
            {

                std::ifstream file;
                file.open(path);

                if (file.is_open())
                {

                    std::string first, line;
                    size_t oid = 0;
                    std::getline(file, first);
                    dataset::BasicArrayObject<O, T> aux = dataset::BasicArrayObject<O, T>();

                    aux.setOID(oid++);

                    std::vector<std::string> vecST = splitString(first, pattern);

                    setDimensionality(std::max(vecST.size(), getDimensionality()));

                    for(size_t x = 0; x < getDimensionality(); x++)
                    {

                        std::string s = vecST[x];
                        T charVec(s.begin(), s.end());
                        aux.set(charVec);

                    }

                    dataset->push_back(aux);

                    while(std::getline(file, line))
                    {

                        dataset::BasicArrayObject<O, T> aux2 = dataset::BasicArrayObject<O, T>();
                        aux2.setOID(oid++);

                        std::regex rgx(pattern);
                        std::sregex_token_iterator iter(line.begin(),line.end(),rgx,-1);
                        std::sregex_token_iterator end;

                        for(; iter != end; ++iter)
                        {

                            std::string st = *iter;
                            T charVec(st.begin(), st.end());
                            aux2.set(charVec);

                        }

                        setDimensionality(std::max(aux2.size(), getDimensionality()));
                        dataset->push_back(aux2);

                    }

                }
                else
                    throw std::runtime_error("Dataset::readText: File not found");

            }
            else
                throw std::runtime_error("Dataset::readText: Invalid parameters");

        }

    public:
        //Constructors and destructor

        Dataset() : dataset(std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>()), dimensionality(0), seed(0) {}

        Dataset(std::vector<dataset::BasicArrayObject<O, T>> _dataset, size_t _dimensionality, size_t _seed = 0, std::string _path = "") :
                dataset(std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(_dataset)), dimensionality(_dimensionality), seed(_seed), path(std::move(_path)) {}

        Dataset(std::string _path, size_t _cardinality, size_t _dimensionality, size_t _seed = 0) : dataset(std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>()), dimensionality(_dimensionality), seed(_seed), path(std::move(_path))
        {
            read(_cardinality, _dimensionality);
        }

        explicit Dataset(std::string _path, std::string pattern="[//s+,;]", size_t _seed = 0) : dataset(std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>()), dimensionality(0), seed(_seed), path(std::move(_path))
        {
            read(0, 0, pattern);
        }

        Dataset(const Dataset& obj) : dataset(std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(*obj.dataset)), dimensionality(obj.dimensionality), seed(obj.seed), path(obj.path) {}

        virtual ~Dataset()
        {

            if (dataset != nullptr)
            {
                for (auto &array : *dataset)
                {
                    array.clear();
                }

                dataset->clear();
                dataset.reset();
            }

        }

        // Virtual methods

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(getSerializedSize());
            size_t offset = 0, cardinality = getCardinality(), serializedSize;

            memcpy(data.get() + offset, &cardinality, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &dimensionality, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, &seed, sizeof(size_t));
            offset += sizeof(size_t);

            size_t pathSize = path.size();
            memcpy(data.get() + offset, &pathSize, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(data.get() + offset, path.c_str(), pathSize);
            offset += pathSize;

            for(size_t x = 0; x < cardinality; x++)
            {

                serializedSize = (*dataset)[x].getSerializedSize();
                memcpy(data.get() + offset, &serializedSize, sizeof(size_t));
                offset += sizeof(size_t);

                memcpy(data.get() + offset, (*dataset)[x].serialize().get(), serializedSize);
                offset += serializedSize;

            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, cardinality, serializedSize;

            memcpy(&cardinality, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&dimensionality, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(&seed, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            size_t pathSize;
            memcpy(&pathSize, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            char* pathChar = new char[pathSize];
            memcpy(pathChar, _data.get() + offset, pathSize);
            offset += pathSize;

            path = std::string(pathChar, pathSize);
            delete[] pathChar;

            if (cardinality == 0)
                return;

            this->clear();
            dataset = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>();

            for(size_t x = 0; x < cardinality; x++)
            {

                memcpy(&serializedSize, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                std::unique_ptr<u_char[]> data = std::make_unique<u_char[]>(serializedSize);
                memcpy(data.get(), _data.get() + offset, serializedSize);
                offset += serializedSize;

                dataset::BasicArrayObject<O, T> aux;
                aux.deserialize(std::move(data));
                dataset->push_back(aux);

            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            size_t size = 0;

            size += sizeof(size_t); // cardinality
            size += sizeof(size_t); // dimensionality
            size += sizeof(size_t); // seed
            size += sizeof(size_t); // path size
            size += path.size(); // path

            for(size_t x = 0; x < this->size(); x++)
            {

                size += sizeof(size_t) + (*dataset)[x].getSerializedSize();

            }

            return size;
        }

        // Setters

        void setDimensionality(size_t _dimensionality)
        {
            dimensionality = _dimensionality;
        }

        void setData(std::vector<dataset::BasicArrayObject<O, T>> _dataset)
        {
            clear();
            dataset = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>(_dataset);
        }

        void setSeed(size_t _seed)
        {
            seed = _seed;
        }

        void setPath(std::string _path)
        {
            path = std::move(_path);
        }

        // Getters

        [[nodiscard]] size_t getDimensionality() const
        {
            return dimensionality;
        }

        std::vector<dataset::BasicArrayObject<O, T>> getData() const
        {
            if (dataset == nullptr)
                return std::vector<dataset::BasicArrayObject<O, T>>();

            return *dataset;
        }

        [[nodiscard]] size_t getSeed() const
        {
            return seed;
        }

        [[nodiscard]] std::string getPath() const
        {
            return path;
        }

        [[nodiscard]] size_t getCardinality() const
        {
            return size();
        }


        // Operators

        friend std::ostream& operator<<(std::ostream& os, const Dataset& obj)
        {

//            os << std::format("Dataset: {0} ; {1} ; {2} ; {3}\n", obj.dataset->size(), obj.dimensionality, obj.seed, obj.path);
//
//            for(size_t x = 0; x < obj.dataset->size(); x++)
//            {
//
//                os << obj.dataset->at(x) << "\n";
//
//            }

            return os;

        }

        bool operator==(const Dataset<O, T>& obj) const
        {

            return isEqual(obj);

        }

        bool operator!=(const Dataset<O, T>& obj) const
        {

            return !isEqual(obj);

        }

        dataset::BasicArrayObject<O, T> operator[](size_t index) const
        {
            if (dataset == nullptr)
                return dataset::BasicArrayObject<O, T>();

            utils::check_range(0, size()-1, index, "Dataset::operator[]: Invalid position");
            return dataset->at(index);
        }

        dataset::BasicArrayObject<O, T>& operator[](size_t index)
        {
            if (dataset == nullptr)
                throw std::runtime_error("Dataset::operator[]: Dataset is empty");

            utils::check_range(0, size()-1, index, "Dataset::operator[]: Invalid position");
            return dataset->at(index);
        }

        // Methods

        void clear()
        {
            if (dataset != nullptr)
            {
                for (auto &array : *dataset)
                {
                    array.clear();
                }

                dataset->clear();
                dataset.reset();
            }
        }

        [[nodiscard]] size_t size() const
        {
            if (dataset == nullptr)
                return 0;

            return dataset->size();
        }

        bool empty()
        {
            return (size() == 0);
        }

        dataset::BasicArrayObject<O, T>& getElement(size_t index)
        {
            if (dataset == nullptr)
                throw std::runtime_error("Dataset::getElement: Dataset is empty");

            utils::check_range(0, size()-1, index, "Dataset::getElement: Invalid position");
            return dataset->at(index);
        }

        std::unique_ptr<Dataset<O, T>> sample(std::variant<size_t, double> size, bool reposition, size_t _seed = 0)
        {

            std::vector<size_t> indexes;
            size_t sz;

            if(size.index() == 0)
                sz = std::get<size_t>(size);
            else
                sz = std::floor(std::get<double>(size) * getCardinality());

            indexes = utils::generateRandomNumbers(0, getCardinality()-1, sz, reposition, (_seed == 0 ? getSeed() : _seed));
            std::vector<dataset::BasicArrayObject<O, T>> aux;

            for(unsigned long & index : indexes)
            {

                aux.push_back((*dataset)[index]);

            }

            return std::make_unique<Dataset<O, T>>(aux, getDimensionality(), (_seed == 0 ? getSeed() : _seed), getPath());

        }

        bool isEqual(const Dataset<O, T>& obj) const
        {

            if(getCardinality() != obj.getCardinality() || getDimensionality() != obj.getDimensionality())
                return false;

            for(size_t x = 0; x < getCardinality(); x++)
            {

                if(operator[](x) != obj[x])
                    return false;

            }

            return true;

        }

        bool contains(dataset::BasicArrayObject<O, T> obj) const
        {

            for(size_t x = 0; x < getCardinality(); x++)
            {

                if(operator[](x) == obj)
                    return true;

            }

            return false;

        }

        void erase(size_t pos)
        {

            if (dataset == nullptr)
                throw std::runtime_error("Dataset::erase: Dataset is empty");

            utils::check_range(0, size()-1, pos, "Dataset::erase: Invalid position");
            dataset->erase(dataset->begin() + pos);

        }

        void erase(dataset::BasicArrayObject<O, T> obj)
        {

            if (dataset == nullptr)
                throw std::runtime_error("Dataset::erase: Dataset is empty");

            for(size_t x = 0; x < getCardinality(); x++)
            {

                if(operator[](x) == obj)
                {

                    erase(x);
                    return;

                }

            }

        }

        void insert(dataset::BasicArrayObject<O, T> obj)
        {

            if (dataset == nullptr)
                dataset = std::make_unique<std::vector<dataset::BasicArrayObject<O, T>>>();

            this->dimensionality = std::max(this->dimensionality, obj.size());
            dataset->push_back(obj);

        }

    };




}


#endif //GERVLIB_DATASET_H
