#!/bin/bash


##########
#
# This macro runs madgraph multiple times, one set for every card found in carddir.
# For each file in there, we replace the Cards/run_card.dat with it and use that name
# as the basename of the output file.
#
##########

carddir="bgruncards"

#how many jobs do you want to runfor each
nJobs=1
#what should the starting index be?  (lets me stack these)
startat=0

#move the original runcard somewhere safe, just in case we break something:
cp Cards/run_card.dat Cards/backup_run_card.dat

for outputname in `ls $carddir/`
  do
      echo "running from $carddir/$outputname"
      cp $carddir/$outputname Cards/run_card.dat
      ./bin/newprocess
      for ((job=0+startat;job<nJobs+startat;job=job+1)); do
	echo " "
	echo "****************************"
	echo "FIRE job  ", $job
	echo "****************************"	
        time ./bin/generate_events 0 $outputname.$job
      done

      #add all the jobs together.  So long as they have the same number of requested events, the weights should combine correctly.
      totJobs=`ls -1 Events/$outputname.*.ttree.root | wc -l`
      hadd Events/sum${nJobs}_$outputname.ttree.root Events/$outputname.*.ttree.root
  done




exit
