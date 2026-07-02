#!/bin/sh
# Amy Campbell
# Run bootstrap analysis after having already run RAxML with however many random starts
source /home/campbela12/miniforge3/bin/activate ST105Env

# alignment input
alignmentfile=$1

# number threads to pass
nthread=$2

# prefix from RAxML run where you just got the best tree from 100 starts
originaltreename=$3

# Name for the tree from the bootstrapping mode run
bootstraptreename=$4

# Name for drawing bipartitions 
partitionname=$5

# number random starts (-# parameter, # bootstraps in this case)
numstarts=$6


raxmlHPC-PTHREADS -T $nthread -m GTRGAMMA -p 19104 -b 19104 -# $numstarts -s $alignmentfile -n $bootstraptreename
#raxmlHPC-PTHREADS -m GTRGAMMA -p 19104 -f b -t "RAxML_bestTree."$originaltreename -z "RAxML_bootstrap."$bootstraptreename -n $partitionname




