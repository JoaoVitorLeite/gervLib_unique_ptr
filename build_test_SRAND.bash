#!/bin/bash

#rm -rf build
#mkdir build
#cd build
#cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SRAND=ON ..
#make
##./gervLib
#ctest --extra-verbose

rm -rf build
mkdir build
# shellcheck disable=SC2164
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_PROCESS_TIME:BOOL=ON -DENABLE_SRAND:BOOL=OFF -DENABLE_TESTS:BOOL=OFF -DCMAKE_C_FLAGS_DEBUG="-g -O3" -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" ..
make

./gervLibTest -INDEX VPTREE -DATASET_TRAIN ../data/cities_norm.csv -DATASET_TRAIN_SEPARATOR , -DATASET_TEST ../data/cities_norm.csv -DATASET_TEST_SEPARATOR , -DISTANCE_FUNCTION EUCLIDEAN -PIVOT_TYPE RANDOM -PIVOT_SAMPLE_SIZE 1.0 -NUM_PIVOTS 2 -SEED 14 -K_MAX 100
