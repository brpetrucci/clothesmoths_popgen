#!/bin/bash

for f in bam_lists/*.txt; do
  donor=`echo $f | cut -d '/' -f 2 | cut -d '_' -f 1`
  angsd -sites allSites_maxhet -b $f  -GL 1 -P 1 $TODO -out $donor
done
