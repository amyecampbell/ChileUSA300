#!/bin/sh
# Amy Campbell
# Make SNIPPY input file from an input folder


outputfile=$1
inputfolder=$2

touch $outputfile

for file in $inputfolder* ; do
	justID=$(basename $file)
	justID=${justID/".fasta"/""}
	printf $justID'\t'$file'\n' >> $outputfile


done
