//
// Created by joaovictor on 21/08/23.
//

#ifndef GERVLIB_NARYTREE_H
#define GERVLIB_NARYTREE_H

#include <vector>
#include <string>
#include <queue>
#include <sstream>
#include <iostream>
#include <memory>

namespace gervLib::naryTree {

    struct NodeNAry {
        std::string value;
        std::vector<std::unique_ptr<NodeNAry>> children;

        NodeNAry() = default;

        explicit NodeNAry(std::string _value) : value(_value) {}
    };

}

#endif //GERVLIB_NARYTREE_H
