#!/bin/bash


##########
#
# This macro runs madgraph multiple times, one set for every card found in carddir.
# For each file in there, we replace the Cards/run_card.dat with it and use that name
# as the basename of the output file.
#
##########

carddir="bgruncards"
nJobs=100

#move the original runcard somewhere safe, just in case we break something:
cp Cards/run_card.dat Cards/backup_run_card.dat

for outputname in `ls $carddir/`
  do
      cp $carddir/$outputname Cards/run_card.dat
      ./bin/newprocess
      for ((job=0;job<nJobs;job=job+1)); do
	echo " "
	echo "****************************"
	echo "FIRE job  ", $iJob
	echo "****************************"	
        time ./bin/generate_events 0 $outputname.$job
      done
      
      hadd Events/sum${nJobs}_$outputname.ttree.root Events/$outputname.*.ttree.root
  done




exit
