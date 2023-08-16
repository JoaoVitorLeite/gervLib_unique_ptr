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

    void check_range(const size_t min, const size_t max, const size_t value, const std::string &message)
    {
        if(min > max)
        {
            throw std::runtime_error("Min value greater than max value | " + message);
        }
        else if(min == max && value != min)
        {
            throw std::runtime_error("Min value equals to max value but value is different from both | " + message);
        }
        else if(value < min || value > max)
        {
            throw std::out_of_range("Value out of range: " + std::to_string(value) + " | " + message);
        }

    }


}

#endif //GERVLIB_UTILS_H
