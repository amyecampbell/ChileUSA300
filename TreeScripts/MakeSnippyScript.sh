#!/bin/sh
# Amy Campbell
# Make script to call to make Snippy alignment (pre-recombination check)
source /home/campbela12/miniforge3/bin/activate ST105Env

refgenome=$1
tabularinput=$2
output=$3

snippy-multi $tabularinput --ref $refgenome  --cpus 16 >  $output


