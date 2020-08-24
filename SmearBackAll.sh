#!/bin/bash

# Convert all .djangoh.txt files through various permutations so they can be smeared and written back as .otree.root files.

dirname=$1
echo dirname: .$dirname.
outdir=$2
echo outdir: .$outdir.
filelist=`ls $dirname/*.djangoh.txt`
mkdir $outdir

degrade=$3

echo taking every .djangoh file in $dirname and processing it with degrade=$degrade into $outdir

for filename in $filelist
do
    outputname=`echo $filename | sed s@$dirname@$outdir@`
    echo running with errors silenced!  root -b -q SmearBackToSimple.C\\\(\\\"$filename\\\",\\\"$outputname\\\",$degrade\\\) 
    root -b -q SmearBackToSimple.C\(\"$filename\",\"$outputname\",$degrade\) 2>/dev/null
done
exit
