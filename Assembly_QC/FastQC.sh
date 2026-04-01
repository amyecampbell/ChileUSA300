#!/bin/sh

source /home/campbela12/miniforge3/bin/activate ST105Env


output=$1
# /scr1/users/campbela12/ST105/FastQC/

input=$2

#/scr1/users/campbela12/ST105/rawreads/ 

#for file in $input*; do
#	fastqc -o $output $file -f fastq

#done

multiqc -n MultiQCST105Raw $output
