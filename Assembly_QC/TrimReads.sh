#!/bin/sh
# Run trim galore on the 99 S. aureus isolates

# e.g., /scr1/users/campbela12/ST105/rawreads/PP.33 to loop through the 
# isolates whose names start with PP.33 (lazy parallelization)
FilePathPrefix=$1

# where to put the trimmed reads
Outputfolder=$2

# fwd suffix (like "1" or "f"); thing that comes before fastq.gz
FwdSuffix=$3

# rev suffix (like "2" or "r"); thing that comes before fastq.gz
RevSuffix=$4

source /home/campbela12/miniforge3/bin/activate ST105Env

for fwdfilename in $FilePathPrefix*$FwdSuffix.fastq.gz ; do


	fwdFileEnd=$FwdSuffix.fastq.gz
	revFileEnd=${fwdFileEnd/$FwdSuffix/$RevSuffix}

	WithUnderscore="_"$fwdFileEnd
	revfilename=${fwdfilename/$fwdFileEnd/$revFileEnd}

	basefilename=$(basename $fwdfilename)
	isolateID=${basefilename/$WithUnderscore/""}

	#  trimming call
	trim_galore --paired --clip_R1 10 --clip_R2 10 --length 50 --basename $isolateID --output_dir  $Outputfolder $fwdfilename $revfilename

done


 

