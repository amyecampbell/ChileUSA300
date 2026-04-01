#!/bin/sh
# Amy Campbell
# Mash search on 99 ST105 isolates

source /home/campbela12/miniforge3/bin/activate ST105Env

# e.g. /scr1/users/campbela12/ST105/contigs/
FilePathPrefix=$1

# outputfolder for the actual output files (e.g. HMM calls) such as /scr1/users/campbela12/ST105/MashOutputs/
outputfolder=$2

# mash db path
mashdb=$3


mkdir -p $outputfolder


for file in $FilePathPrefix* ; do 

	inputname=$(basename $file)
	output="_Mash.tab"
	inputextension=".fasta"
	outputfilename=${inputname/$inputextension/$output}
	outputpath=$outputfolder$outputfilename
	mash screen -w -p 4 $mashdb $file  > $outputpath
	
	## commenting this out for now
	#cleaned=${outputpath/".tab"/"_top.tab"}
	#sort -grk2 $outputpath | head > $cleaned
	# rm $outputpath
done
