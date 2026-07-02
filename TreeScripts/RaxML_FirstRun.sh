#!/bin/sh
# Amy Campbell
#  make ML tree from Snippy alignment
source /home/campbela12/miniforge3/bin/activate ST105Env

# alignment input
alignmentfile=$1

# number threads to pass
nthread=$2

# prefix for raxml
name=$3

# number random starts (-# parameter)
numstarts=$4

raxmlHPC-PTHREADS -T $nthread -m GTRGAMMA -p 19104 -# $numstarts -s $alignmentfile -n $name




