//
// Created by joaovictor on 15/08/23.
//

#ifndef GERVLIB_CONFIGURE_H
#define GERVLIB_CONFIGURE_H

#include "Utils.h"
#include <iostream>
#include <variant>
#include <gmpxx.h>
#include <type_traits>

namespace gervLib::configure {

    std::string baseOutputPath = SOURCE_OUTPUT_PATH;

    unsigned long long IORead = 0, IOWrite = 0;

    size_t oid_size = 8, double_size = 20, str_size = 25;

    template <typename... Ti>
    std::string variant2string(std::variant<Ti...> var)
    {
        return std::visit([]<typename T>(const T& e){
            if constexpr (std::is_same_v<T, std::string>) {
                return e;
            } else { // float/int
                return std::to_string(e);
            }
        }, var);
    }

    void configure()
    {

#ifdef ENABLE_DEBUG
        std::cout << "\nConfiguring global variables...\n\n" << std::endl;
#endif

        gervLib::utils::createFolderIfNotExists(baseOutputPath);

    }


}

#endif //GERVLIB_CONFIGURE_H
