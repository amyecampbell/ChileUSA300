#!/bin/sh
# Amy Campbell
# Sequence type isolates using mlst package (T. Seemann)

source /home/campbela12/miniforge3/bin/activate ST105Env

# path with all the contigs (e.g.)  /scr1/users/campbela12/ST105/ST105Tree_11_2023/SequenceFiles/ 
FilePathPrefix=$1

# where to put the output (e.g. /scr1/users/campbela12/ST105/ST105Tree_11_2023/)
Outputfolder=$2

# Name for the output file (e.g., ST105treeMLST.txt
OutputName=$3 

ulimit -s 32768

mlst --scheme saureus $FilePathPrefix* > $Outputfolder$OutputName

