#!/bin/sh
# Amy Campbell
# Recombination  w clonalframe ML
source /home/campbela12/miniforge3/bin/activate ST105Env

Treefile=$1
Alignment=$2
Prefix=$3
#OutputAlign=$4

ClonalFrameML $Treefile $Alignment $Prefix -output_filtered true -emsim 100

# not doing this anymore since adding -output_filtered true
#/home/campbela12/Documents/software/cfml-maskrc/cfml-maskrc.py --aln $Alignment --out $OutputAlign --regions $Prefix"_regions" $Prefix
