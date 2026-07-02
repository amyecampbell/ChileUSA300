#!/bin/sh
# Amy Campbell
# Make SNIPPY input file for the full tree of ST105 isolates


outputfile=$1
inputfolder=$2

touch $outputfile

for file in $inputfolder* ; do
	justID=$(basename $file)
	justID=${justID/".fasta"/""}
	printf $justID'\t'$file'\n' >> $outputfile


done
