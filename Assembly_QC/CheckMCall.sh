#!/bin/sh
# Amy Campbell
# CheckM QC on 99 ST105 isolates

source /home/campbela12/miniforge3/bin/activate ST105Env

# e.g. /scr1/users/campbela12/ST105/contigs/
FilePathPrefix=$1

# outputfolder for the actual output files (e.g. HMM calls) such as /scr1/users/campbela12/ST105/CheckMoutput/
outputfolder=$2

# marker file /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers
markerfile=$3

summaryoutput=$4


mkdir -p $outputfolder

#checkm taxon_set species  "Staphylococcus aureus" /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers

checkm analyze -t 16 -x fasta $markerfile $FilePathPrefix $outputfolder

checkm qa $markerfile $outputfolder > $summaryoutput


