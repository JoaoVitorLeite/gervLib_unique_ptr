#!/bin/bash

outputPath="output_files"
index="VPTREE" ####
datasetName="CITIES" ####
datasetTrain="../data/cities_norm.csv" ####
datasetTrainSeparator=","
datasetTest="../data/cities_norm.csv" ####
datasetTestSeparator=","
distanceFunction="EUCLIDEAN"
pivotType=("BPP" "CONVEX" "FFT" "SSS" "MAXSEPARATED" "MAXVARIANCE" "PCA" "IS" "HFI" "WDR" "SELECTION" "KMEDOIDS" "RANDOM")
sampleSize=("1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0") ####
numPivots="2" ####
# shellcheck disable=SC2207
seed=($(shuf -i 0-500000 -n 13))
kMax="100"
numPerLeaf=55 ####
numBins="256"
pageSize="0"
numQueryPerFile="10000"
storePivotLeaf="TRUE"
storeDirectoryNode="FALSE"
storeLeafNode="TRUE"
useLAESA="TRUE" ####

sed -i '104s/.*/static const int leafslots = '${numPerLeaf}';/' libs/bptree/btree.h

rm -rf build
mkdir build
# shellcheck disable=SC2164
cd build
mkdir $outputPath
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_PROCESS_TIME:BOOL=ON -DENABLE_SRAND:BOOL=OFF -DENABLE_TESTS:BOOL=OFF -DCMAKE_C_FLAGS_DEBUG="-g -O3" -DCMAKE_CXX_FLAGS_DEBUG="-g -O3" .. -DSOURCE_OUTPUT_PATH:STRING=$outputPath
make

for((i=0; i<2; i++));
do
    nohup ./gervLibTest -INDEX ${index} -DATASET_NAME ${datasetName} -DATASET_TRAIN ${datasetTrain} -DATASET_TRAIN_SEPARATOR ${datasetTrainSeparator} -DATASET_TEST ${datasetTest} -DATASET_TEST_SEPARATOR ${datasetTestSeparator} -DISTANCE_FUNCTION ${distanceFunction} -PIVOT_TYPE "${pivotType[$i]}" -PIVOT_SAMPLE_SIZE "${sampleSize[$i]}" -NUM_PIVOTS ${numPivots} -SEED "${seed[$i]}" -K_MAX ${kMax} -NUM_PER_LEAF ${numPerLeaf} -NUM_BINS ${numBins} -PAGE_SIZE ${pageSize} -NUM_QUERIES_PER_FILE ${numQueryPerFile} -STORE_PIVOT_LEAF ${storePivotLeaf} -STORE_DIRECTORY_NODE ${storeDirectoryNode} -STORE_LEAF_NODE ${storeLeafNode} -USE_LAESA ${useLAESA} &
done
