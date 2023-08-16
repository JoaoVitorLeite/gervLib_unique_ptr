//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_DATASETWRAPPER_H
#define GERVLIB_DATASETWRAPPER_H

#include <map>
#include "Dataset.h"

namespace gervLib::dataset
{

    template <typename O, typename T>
    class DatasetWrapper
    {

    private:
        std::map<O, O> index;

    private:
        void mapIndex(std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {

            for(size_t i = 0; i < dataset->getCardinality(); i++)
            {

                index.insert(std::make_pair(static_cast<O>(i), dataset->getElement(i).getOID()));
                dataset->getElement(i).setOID(static_cast<O>(i));

            }

        }

    public:
        DatasetWrapper() = default;

        explicit DatasetWrapper(std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {
            mapIndex(dataset);
        }

        ~DatasetWrapper()
        {

            index.clear();

        }

        void resetIndex(std::unique_ptr<dataset::Dataset<O, T>>& dataset)
        {

            for(size_t i = 0; i < dataset->getCardinality(); i++)
            {

                dataset->getElement(i).setOID(index[static_cast<O>(i)]);

            }

        }

        void resetIndex(dataset::BasicArrayObject<O, T>& element)
        {

            element.setOID(index[element.getOID()]);

        }


    };


}

#endif //GERVLIB_DATASETWRAPPER_H
