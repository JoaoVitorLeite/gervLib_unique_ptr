//
// Created by joaovictor on 15/08/23.
//

#ifndef GERVLIB_UTILS_H
#define GERVLIB_UTILS_H

#include <string>
#include <filesystem>

namespace gervLib::utils
{

    void createFolderIfNotExists(const std::string &path)
    {

        if (!std::filesystem::exists(path)) {
            std::filesystem::create_directory(path);
        }

    }

}

#endif //GERVLIB_UTILS_H
