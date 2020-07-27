#!/bin/bash

# Convert all .djangoh.txt files through various permutations so they can be smeared and written back as .otree.root files.

dirname=${$1%/}
echo dirname: .$dirname.
filelist=`ls ./$dirname/*.djangoh.txt`


for filename in $filelist
  do
      echo running with errors silenced!  root -b -q SmearBackToSimple.C\\\(\\\"$filename\\\"\\\) 
      root -b -q SmearBackToSimple.C\(\"$filename\"\) 2>/dev/null
done
exit
