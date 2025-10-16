#!/bin/bash

for f1 in angsd_out/*; do
  donor1=`echo $f1 | cut -d '/' -f 2 | cut -d '.' -f 1 | tr -cd '0-9'`
  
  realSFS $f1 >donor${donor1}.sfs

  for f2 in angsd_out/*; do
    # Skip if the files are the same
    [ "$f1" = "$f2" ] && continue
    # Skip symmetric duplicates by comparing names
    [[ "$f1" > "$f2" ]] && continue
    
    donor2=`echo $f2 | cut -d '/' -f 2 | cut -d '.' -f 1 | tr -cd '0-9'`

    filename="p${donor1}_${donor2}"

    realSFS $f1 $f2 -P 24 > $filename.sfs
    
    realSFS fst index $f1 $f2 -sfs $filename.sfs -fstout $filename
    
    realSFS fst stats $filename.fst.idx > $filename.fst.stats
    
    done
done

