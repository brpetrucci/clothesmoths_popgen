#!/bin/bash

datasets=("12p" "50p" "mass")
numdonors=(25 23 7)

for i in "${!datasets[@]}"; do
  for K in `seq 2 ${numdonors[$i]}`; do
      for rep in `seq 1 100`; do
        ./NGSadmix -likes data_raw/${datasets[$i]}.beagle.gz -K $K -P 10 -o ${datasets[$i]}/ngsadmix_k${K}.${rep};
      done
  done
done
