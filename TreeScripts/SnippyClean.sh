#!/bin/sh
# Amy Campbell
# Make script to call to clean up Snippy alignment (pre-recombination check)

source /home/campbela12/miniforge3/bin/activate ST105Env

Alignment=$1
Output=$2

snippy-clean_full_aln $Alignment > $Output



