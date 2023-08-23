//
// Created by joaovictor on 21/08/23.
//

#ifndef GERVLIB_INDEXFACTORY_H
#define GERVLIB_INDEXFACTORY_H

#include "SequentialScan.h"
#include "LAESA.h"

namespace gervLib::index
{

    template <typename O, typename T>
    class IndexFactory
    {
    public:
        static std::unique_ptr<Index<O, T>> createIndex(INDEX_TYPE type)
        {
            if (type == SEQUENTIAL_SCAN_t)
                return std::make_unique<SequentialScan<O, T>>();
            else if (type == LAESA_t)
                return std::make_unique<LAESA<O, T>>();
            else
                throw std::invalid_argument("Invalid INDEX_TYPE");
        }
    };

}

#endif //GERVLIB_INDEXFACTORY_H
