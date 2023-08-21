cmake_minimum_required(VERSION 3.25)
project(gervLib
        VERSION 1.0
        LANGUAGES CXX
        DESCRIPTION "Similarity search Library")

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_COMPILER gcc)

# Set build options
set(ENABLE_TESTS OFF CACHE BOOL "Enable tests")

# Timer
set(ENABLE_PROCESS_TIME OFF CACHE BOOL "Enable timer")

# Set output path for generated source files
set(SOURCE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/output_dir" CACHE STRING "Output path for generated source files")
add_compile_definitions(SOURCE_OUTPUT_PATH="${SOURCE_OUTPUT_PATH}")

# Set build options
set(ENABLE_SRAND OFF CACHE BOOL "Enable srand() for random number generation")

# Set build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

# Set build type options
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DENABLE_DEBUG)
    add_definitions(-DENABLE_PROGRESS_BAR)
    set(ENABLE_TESTS ON)
endif ()

# Set test options
if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

# Set build type options
if (ENABLE_SRAND)
    add_definitions(-DENABLE_SRAND)
endif ()

if (ENABLE_PROCESS_TIME)
    add_definitions(-DENABLE_PROCESS_TIME)
endif ()

# Add include directories
include_directories(src)
include_directories(libs/uuid)
include_directories(libs/eigen)
include_directories(include/configure)
include_directories(include/utils)
include_directories(include/dataset)
include_directories(include/serialize)
include_directories(include/distance)
include_directories(include/factory)
include_directories(include/pivots)
include_directories(include/kmedoids)
include_directories(include/indexes)
include_directories(include/memory)
include_directories(include/query)
include_directories(include/matrix)

# Add source files
add_executable(gervLib
        src/main.cpp
        libs/uuid/libuuidpp.hpp
        include/configure/Configure.h
        include/utils/Utils.h
        include/dataset/BasicArrayObject.h
        include/serialize/Serialize.h
        include/dataset/Dataset.h
        include/dataset/DatasetWrapper.h
        include/distance/DistanceFunction.h
        include/distance/EuclideanDistance.h
        include/distance/EditDistance.h
        include/factory/DistanceFactory.h
        include/pivots/Pivot.h
        include/pivots/RandomPivots.h
        include/pivots/FFTPivots.h
        include/pivots/ConvexPivots.h
        include/pivots/MaxSeparatedPivots.h
        include/pivots/MaxVariancePivots.h
        include/pivots/SelectionPivots.h
        include/pivots/PCAPivots.h
        include/pivots/SSSPivots.h
        include/pivots/HFIPivots.h
        include/pivots/ISPivots.h
        include/pivots/WDRPivots.h
        include/pivots/BPPPivots.h
        include/kmedoids/Cluster.h
        include/kmedoids/Kmedoids.h
        include/pivots/KmedoidsPivots.h
        include/factory/PivotFactory.h
        include/factory/PivotFactory.h
        include/indexes/Index.h
        include/memory/Page.h
        include/memory/PageManager.h
        include/query/Query.h
        include/indexes/SequentialScan.h
        include/matrix/Matrix.h
        include/indexes/LAESA.h
)

# Add libraries
if (WIN32)
    target_link_libraries (gervLib libuuid.dll)
else()
    target_link_libraries (gervLib libuuid.so)
endif()

# Add compile options
add_compile_options("-std=c++20")

# Enable C++20 features
target_compile_features(gervLib PRIVATE cxx_std_20)

# Set common compiler flags
target_compile_options(gervLib PRIVATE
        -Wall
        -Wextra
        -pedantic
        -pipe
        -g
        -std=gnu++2a
        -fPIC
        -D_REENTRANT
        # Add more flags as needed
)