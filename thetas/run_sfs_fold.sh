#!/bin/bash

for file in donor*.saf.idx; do
  N=$(echo "$file" | sed -E 's/donor([0-9]+)\.saf\.idx/\1/')

  realSFS donor${N}.saf.idx -fold 1 > donor${N}.sfs

  realSFS saf2theta donor${N}.saf.idx -sfs donor${N}.sfs -outname donor${N} -fold 1

  thetaStat print donor${N}.thetas.idx 2>/dev/null |head > donor${N}.thetas.print

  thetaStat do_stat donor${N}.thetas.idx 

  thetaStat do_stat donor${N}.thetas.idx -win 50000 -step 10000 -outnames donor${N}.thetasWindow.gz
done
