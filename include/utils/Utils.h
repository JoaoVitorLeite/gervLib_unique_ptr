//
// Created by joaovictor on 15/08/23.
//

#ifndef GERVLIB_UTILS_H
#define GERVLIB_UTILS_H

#include <string>
#include <filesystem>
#include <vector>
#include <set>
#include <random>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <sys/time.h>
#include <sys/resource.h>
#include <fstream>
#include <iostream>
#include "libuuidpp.hpp"

namespace gervLib::utils
{

    std::string generateExperimentID()
    {
        return libuuidpp::uuid::random().string();
    }

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

    template<typename generator_type = size_t>
    struct Random {

        std::mt19937_64 generator;

        explicit Random(generator_type seed = 0) {
#ifndef ENABLE_SRAND
            generator.seed(seed);
#else
            srand(seed);
#endif

        }

        generator_type operator()(generator_type min, generator_type max) {
#ifndef ENABLE_SRAND
            std::uniform_int_distribution<generator_type> distribution(min, max);
            return distribution(generator);
#else
            return min + (rand() % (max - min + 1));
#endif
        }

        generator_type operator()(generator_type max) {
            return (*this)(0, max);
        }

        generator_type operator()() {
#ifndef ENABLE_SRAND
            return (*this)(0, std::numeric_limits<generator_type>::max());
#else
            return (*this)(0, RAND_MAX);
#endif
        }

    };

    std::vector<size_t> generateRandomNumbers(size_t min, size_t max, size_t size, bool reposition, size_t seed)
    {

        std::vector<size_t> numbers;
        Random<size_t> random(seed);

        if(reposition)
        {

            for(size_t i = 0; i < size; i++)
            {
                numbers.push_back(random(min, max));
            }

        }
        else
        {

            if(max - min + 1 < size)
            {
                throw std::runtime_error("The range of numbers is smaller than the size of the vector");
            }

            std::set<size_t> aux;

            while(aux.size() < size)
            {
                aux.insert(random(min, max));
            }

            for(auto &i : aux)
            {
                numbers.push_back(i);
            }

        }

        return numbers;

    }

    struct Timer
    {

#ifdef ENABLE_PROCESS_TIME
        struct timeval start_time_u, start_time_s, end_time_u, end_time_s;
        struct rusage usage;
#else
        std::chrono::high_resolution_clock::time_point start_time;
        std::chrono::high_resolution_clock::time_point end_time;
#endif

    public:
        void start()
        {
#ifdef ENABLE_PROCESS_TIME
            getrusage(RUSAGE_SELF, &usage);
            start_time_u = usage.ru_utime;
            start_time_s = usage.ru_stime;
#else
            start_time = std::chrono::high_resolution_clock::now();
#endif
        }

        void stop()
        {
#ifdef ENABLE_PROCESS_TIME
            getrusage(RUSAGE_SELF, &usage);
            end_time_u = usage.ru_utime;
            end_time_s = usage.ru_stime;
#else
            end_time = std::chrono::high_resolution_clock::now();
#endif
        }

        [[nodiscard]] long long getElapsedTime() const
        {
#ifdef ENABLE_PROCESS_TIME
            return (end_time_u.tv_sec - start_time_u.tv_sec) * 1000000L + (end_time_u.tv_usec - start_time_u.tv_usec) +
                   (end_time_s.tv_sec - start_time_s.tv_sec) * 1000000L + (end_time_s.tv_usec - start_time_s.tv_usec);
#else
            return std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
#endif
        }

        [[nodiscard]] long long getElapsedTimeUser() const
        {
#ifdef ENABLE_PROCESS_TIME
            return (end_time_u.tv_sec - start_time_u.tv_sec) * 1000000L + (end_time_u.tv_usec - start_time_u.tv_usec);
#else
            return -1;
#endif
        }

        [[nodiscard]] long long getElapsedTimeSystem() const
        {
#ifdef ENABLE_PROCESS_TIME
            return (end_time_s.tv_sec - start_time_s.tv_sec) * 1000000L + (end_time_s.tv_usec - start_time_s.tv_usec);
#else
            return -1;
#endif
        }

    };

    std::vector<std::string> searchPathContainingSubstring(const std::string& directoryPath, const std::string& substring) {

        std::vector<std::string> ans;

        for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
            if (std::filesystem::is_directory(entry) && entry.path().filename().string().find(substring) != std::string::npos) {
                ans.push_back(entry.path().string());
            }
        }

        return ans;
    }

    std::vector<std::string> searchFilesWithSubstring(const std::string& path, const std::string& substring)
    {

        std::vector<std::string> ans;

        for(const auto & entry : std::filesystem::directory_iterator(path))
        {
            if(std::filesystem::is_regular_file(entry) && entry.path().filename().string().find(substring) != std::string::npos)
            {
                ans.push_back(entry.path().string());
            }
        }

        return ans;

    }

    std::string generatePathByPrefix(const std::string& path, const std::string& prefix, bool createDir = true)
    {

        std::vector<std::string> paths = searchPathContainingSubstring(path, prefix);

        if(paths.empty())
        {

            std::filesystem::path p(path);
            p /= prefix;
            p += "_0";

            if(createDir)
            {

                if(!std::filesystem::create_directory(p))
                {
                    throw std::runtime_error("generatePathByPrefix(): Could not create directory: " + p.string());
                }

            }

            return p.string();

        }
        else
        {

            std::sort(paths.begin(), paths.end(), [](const std::string& path1, const std::string& path2) -> bool
            {

                std::string dirName1 = std::filesystem::path(path1).filename().string();
                std::string dirName2 = std::filesystem::path(path2).filename().string();

                size_t pos1 = dirName1.find_last_not_of("0123456789");
                size_t pos2 = dirName2.find_last_not_of("0123456789");

                if (pos1 == std::string::npos) {
                    return false;
                }

                if (pos2 == std::string::npos) {
                    return true;
                }

                int lastNumber1 = std::stoi(dirName1.substr(pos1 + 1));
                int lastNumber2 = std::stoi(dirName2.substr(pos2 + 1));

                return lastNumber1 > lastNumber2;

            });

            std::filesystem::path p(paths[0]);
            std::string dirName = p.filename().string();
            size_t pos = dirName.find_last_not_of("0123456789");
            int lastNumber = std::stoi(dirName.substr(pos + 1));
            p.replace_filename(prefix + "_" + std::to_string(lastNumber + 1));

            if(createDir)
            {

                if(!std::filesystem::create_directory(p))
                {
                    throw std::runtime_error("generatePathByPrefix(): Could not create directory: " + p.string());
                }

            }

            return p.string();

        }

    }

    void deleteDirectory(const std::filesystem::path& path) {
        try {
            if (std::filesystem::exists(path) && std::filesystem::is_directory(path)) {
                for (const auto& entry : std::filesystem::directory_iterator(path)) {
                    if (std::filesystem::is_directory(entry)) {
                        deleteDirectory(entry.path());
                    } else {
                        std::filesystem::remove(entry.path());
                    }
                }
                std::filesystem::remove(path); // Remove the directory itself
            }
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    std::string generateFileByPrefix(const std::string& path, const std::string& prefix, bool createFile = true, const std::string& extension = ".txt")
    {

        std::vector<std::string> files = searchFilesWithSubstring(path, prefix);

        if(files.empty())
        {
            std::filesystem::path p(path);
            std::string fileName = (p / (prefix + "_0" + extension)).string();

            if(createFile)
            {
                std::ofstream file(fileName);

                if(!file.is_open())
                {
                    throw std::runtime_error("generateFileByPrefix(): Could not create file: " + fileName);
                }

                file.close();
            }

            return fileName;
        }
        else
        {

            std::sort(files.begin(), files.end(), [](const std::string& filePath1, const std::string& filePath2) -> bool
            {

                std::string fileName1 = std::filesystem::path(filePath1).filename().stem().string();
                std::string fileName2 = std::filesystem::path(filePath2).filename().stem().string();

                size_t pos1 = fileName1.find_last_not_of("0123456789");
                size_t pos2 = fileName2.find_last_not_of("0123456789");

                if (pos1 == std::string::npos) {
                    return false;
                }

                if (pos2 == std::string::npos) {
                    return true;
                }

                int lastNumber1 = std::stoi(fileName1.substr(pos1 + 1));
                int lastNumber2 = std::stoi(fileName2.substr(pos2 + 1));

                return lastNumber1 > lastNumber2;

            });

            std::string fileName = std::filesystem::path(files[0]).filename().stem().string();
            size_t pos = fileName.find_last_not_of("0123456789");
            int lastNumber = std::stoi(fileName.substr(pos + 1));

            std::filesystem::path p(path);
            std::string gen_fileName = (p / (prefix + "_" + std::to_string(lastNumber + 1) + extension)).string();

            if(createFile)
            {
                std::ofstream file(gen_fileName);

                if(!file.is_open())
                {
                    throw std::runtime_error("generateFileByPrefix(): Could not create file: " + gen_fileName);
                }

                file.close();
            }

            return gen_fileName;

        }

    }

    std::string getPathByPrefix(const std::string& path, const std::string& prefix)
    {

        std::vector<std::string> paths = searchPathContainingSubstring(path, prefix);

        if(paths.empty())
        {

            throw std::runtime_error("getPathByPrefix(): Could not find path with prefix: " + prefix);

        }
        else
        {

            std::sort(paths.begin(), paths.end(), [](const std::string& path1, const std::string& path2) -> bool
            {

                std::string dirName1 = std::filesystem::path(path1).filename().string();
                std::string dirName2 = std::filesystem::path(path2).filename().string();

                size_t pos1 = dirName1.find_last_not_of("0123456789");
                size_t pos2 = dirName2.find_last_not_of("0123456789");

                if (pos1 == std::string::npos) {
                    return false;
                }

                if (pos2 == std::string::npos) {
                    return true;
                }

                int lastNumber1 = std::stoi(dirName1.substr(pos1 + 1));
                int lastNumber2 = std::stoi(dirName2.substr(pos2 + 1));

                return lastNumber1 > lastNumber2;

            });

            std::filesystem::path p(paths[0]);

            return p.string();

        }

    }

    std::string getFileByPrefix(const std::string& path, const std::string& prefix)
    {

        std::vector<std::string> files = searchFilesWithSubstring(path, prefix);

        if(files.empty())
        {
            throw std::runtime_error("getFileByPrefix(): Could not find file with prefix: " + prefix);
        }
        else
        {

            std::sort(files.begin(), files.end(), [](const std::string& filePath1, const std::string& filePath2) -> bool
            {

                std::string fileName1 = std::filesystem::path(filePath1).filename().stem().string();
                std::string fileName2 = std::filesystem::path(filePath2).filename().stem().string();

                size_t pos1 = fileName1.find_last_not_of("0123456789");
                size_t pos2 = fileName2.find_last_not_of("0123456789");

                if (pos1 == std::string::npos) {
                    return false;
                }

                if (pos2 == std::string::npos) {
                    return true;
                }

                int lastNumber1 = std::stoi(fileName1.substr(pos1 + 1));
                int lastNumber2 = std::stoi(fileName2.substr(pos2 + 1));

                return lastNumber1 > lastNumber2;

            });

            std::string fileName = std::filesystem::path(files[0]).filename().string();
            return fileName;

        }

    }


}

#endif //GERVLIB_UTILS_H
