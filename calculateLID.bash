#!/bin/bash

outputPath="output_files"

rm -rf build
mkdir build
# shellcheck disable=SC2164
cd build
mkdir $outputPath
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_PROCESS_TIME:BOOL=ON -DENABLE_SRAND:BOOL=OFF -DENABLE_TESTS:BOOL=OFF -DCMAKE_C_FLAGS_DEBUG="-g -O3" -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" .. -DSOURCE_OUTPUT_PATH:STRING=$outputPath
make

num_procs=$(nproc --all)
num_procs=$((num_procs-9)) ####
split --number=l/$num_procs -d ../data/cities.csv ../data/cities_part --additional-suffix=.csv ####
filePath=($(ls ../data/cities_part*.csv)) ####
fileName=($(ls ../data/cities_part*.csv | xargs -n 1 basename)) ####
len=${#filePath[@]}

./gervLibSerialize -DATASET ${filePath[$i]} -DATASET_SEPARATOR , -PIVOT_SAMPLE_SIZE 0.5 -NUM_PIVOTS 2 -NUM_PER_LEAF 55 -FILE_NAME ${fileName[$i]} -OUTPUT_PATH ../data

for ((i=0; i<$len; i++))
do
    ./gervLibLID -INDEX_FOLDER output_files/VPTREE_MAX_VAR_SERIALIZE -DATASET ${filePath[$i]} -DATASET_SEPARATOR , -FILE_NAME ${fileName[$i]} -OUTPUT_PATH ../data
done

cat ../data/cities_part*lid* > ../data/cities_lid.csv ####
python3 ../scripts/LidSplit.py ../data/cities_lid.csv ../data/cities.csv ####
