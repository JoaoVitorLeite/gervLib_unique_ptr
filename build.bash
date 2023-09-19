#!/bin/bash

outputPath="output_files"

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_PROCESS_TIME:BOOL=ON -DENABLE_SRAND:BOOL=OFF -DENABLE_TESTS:BOOL=OFF -DCMAKE_C_FLAGS_DEBUG="-g -O3" -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" .. -DSOURCE_OUTPUT_PATH:STRING=$outputPath
make

