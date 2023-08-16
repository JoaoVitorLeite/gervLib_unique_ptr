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

}

#endif //GERVLIB_UTILS_H
