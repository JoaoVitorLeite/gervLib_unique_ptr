//
// Created by joaovictor on 17/08/23.
//

#ifndef GERVLIB_INDEX_H
#define GERVLIB_INDEX_H

#include "Dataset.h"
#include "DistanceFactory.h"
#include "PivotFactory.h"

namespace gervLib::index
{

    enum INDEX_TYPE {SEQUENTIAL_SCAN, LAESA, VPTREE, UNKNOWN};

    enum MEMORY_STATUS {IN_MEMORY, IN_DISK, NONE};

    std::map<INDEX_TYPE, std::string> indexTypeMap = {
            {SEQUENTIAL_SCAN, "SEQUENTIAL_SCAN"},
            {LAESA, "LAESA"},
            {VPTREE, "VPTREE"}
    };

    std::map<std::string, INDEX_TYPE> indexTypeMapReverse = {
            {"SEQUENTIAL_SCAN", SEQUENTIAL_SCAN},
            {"LAESA", LAESA},
            {"VPTREE", VPTREE}
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
    class Index
    {



    };


}

#endif //GERVLIB_INDEX_H
