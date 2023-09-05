//
// Created by joaoleite on 8/30/23.
//

#ifndef GERVLIB_EQUIDEPTH_H
#define GERVLIB_EQUIDEPTH_H

#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <gmpxx.h>
#include <type_traits>
#include "Configure.h"

namespace gervLib::equidepth
{

    class Bins
    {
    public:
        long id;
        double min, max;

        Bins()
        {

            id = 0;
            min = std::numeric_limits<double>::max();
            max = std::numeric_limits<double>::lowest();

        }

        Bins(long id_)
        {

            id = id_;
            min = std::numeric_limits<double>::max();
            max = std::numeric_limits<double>::lowest();

        }

        Bins(long id_, double min_, double max_)
        {

            id = id_;
            min = min_;
            max = max_;

        }

        ~Bins(){

        }


        bool isInterval(double test)
        {

            return ((test >= min) && (test <= max));

        }

    };

    template <class T>
    class EquiDepth
    {
    private:
        long num_bins;
        size_t pivot_num;
        bool toINF;
        std::vector<std::vector<Bins>> bins;
        bool inMemory;

    public:
        EquiDepth()
        {

            num_bins = -1;
            pivot_num = -1;
            toINF = true;
            bins = std::vector<std::vector<Bins>>();
            inMemory = true;

        }

        EquiDepth(long num_bins_, size_t pivot_num_, bool toINF_ = true)
        {

            num_bins = num_bins_;
            pivot_num = pivot_num_;
            toINF = toINF_;
            bins = std::vector<std::vector<Bins>>(pivot_num, std::vector<Bins>(num_bins));
            inMemory = true;

        }

        void clear()
        {

            for(auto v : bins)
            {

                v.clear();

            }
            bins.clear();

            inMemory = false;

        }

        ~EquiDepth()
        {

            clear();

        }

        void build(std::vector<std::vector<double>>& aux)
        {

            std::vector<double> v(aux.size());

            for(size_t i = 0; i < pivot_num; i++)
            {

                for(size_t k = 0; k < aux.size(); k++)
                    v[k] = aux[k][i];

                std::sort(v.begin(), v.end());
                size_t index = 0;
                size_t size = v.size()/num_bins;
                size_t sz = 0;
                size_t j = 0;

                //cout << size << endl;
                double last_value = -1.0;

                for(j = 0; j < v.size(); j++)
                {

                    if(sz >= size && v[j] != last_value)
                    {

                        if(index != (num_bins - 1))
                        {

                            bins[i][index+1].min = bins[i][index].max;
                            bins[i][index].max = std::nextafter(bins[i][index].max, std::numeric_limits<double>::lowest());
                            index++;

                        }

                        sz = 0;

                    }

                    bins[i][index].min = std::min(bins[i][index].min, v[j]);
                    bins[i][index].max = std::max(bins[i][index].max, v[j]);
                    sz++;
                    last_value = v[j];

                }

                if(toINF)
                {

                    bins[i][num_bins-1].max = std::numeric_limits<double>::max();

                }


            }

        }

        void print()
        {

            std::cout << "PIVOT NUM: " << bins.size() << std::endl;
            std::cout << "NUM BINS: " << bins[0].size() << std::endl << std::endl;

            for(size_t i = 0; i < bins.size(); i++)
            {

                for(size_t j = 0; j < bins[0].size(); j++)
                {

                    std::cout << "[" << bins[i][j].min << ", " << bins[i][j].max << ")" << std::endl;

                }

                std::cout << "\n\n";

            }

        }

        void saveToFile(const std::string& folder)
        {

            std::ofstream file(folder + std::filesystem::path::preferred_separator + "equi_depth.txt");
            file << bins.size() << " "  << bins[0].size() << std::endl;

            for(size_t i = 0; i < bins.size(); i++)
            {

                for(size_t j = 0; j < bins[0].size(); j++)
                {

                    file << bins[i][j].min << " " << bins[i][j].max << std::endl;

                }

            }

            file.close();

        }

        void readFromFile(const std::string& folder)
        {

            std::ifstream file(folder + std::filesystem::path::preferred_separator + "equi_depth.txt");
//            size_t pivot_num;
            double minV, maxV;

            file >> pivot_num;
            file >> num_bins;

            bins.resize(pivot_num, std::vector<Bins>(num_bins));

            for(size_t i = 0; i < pivot_num; i++)
            {

                for(long j = 0; j < num_bins; j++)
                {

                    file >> minV;
                    file >> maxV;

                    bins[i][j] = Bins(j, minV, maxV);

                }

            }

            file.close();

        }

        long getNumberOfBins()
        {

            return num_bins;

        }

        long getBin(size_t pivot, double value)
        {

            for(long i = 0; i < num_bins; i++)
            {

                if(bins[pivot][i].isInterval(value))
                {

                    return i;

                }

            }

            throw std::runtime_error("Bins definition error !_!");

        }

        std::pair<double, double> getInterval(size_t pivot, long id)
        {

            return std::make_pair(bins[pivot][id].min, bins[pivot][id].max);

        }

        void load(const std::string& folder)
        {

            if(!inMemory)
            {

                readFromFile(folder);
                inMemory = true;

            }

        }

    };

}

#endif //GERVLIB_EQUIDEPTH_H
