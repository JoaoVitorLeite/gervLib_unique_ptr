#!/bin/bash

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SRAND=ON ..
make
#./gervLib
ctest --extra-verbose