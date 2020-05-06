#!/bin/bash

# Convert all .djangoh.txt files through various permutations so they can be smeared and written back as .otree.root files.

dirname=${$1%/}
echo dirname: .$dirname.
filelist=`ls ./$dirname/*.djangoh.txt`


for filename in $filelist
  do
      root -b -q SmearBackToSimple.C\(\"$filename\"\)
done
exit
