#!/bin/sh
# Amy Campbell
# Assemble ST105 Genomes from trimmed reads

source /home/campbela12/miniforge3/bin/activate ST105Env

# path that contains the input reads in name_val_1.fq.gz (for fwd) 
# form. E.g., /scr1/users/campbela12/ST105/TrimmedReads/PP.33 for
# isolates whose names start with PP.33 (lazy parallelization)

FilePathPrefix=$1

# where to put the shovill output, generally
Outputfolder=$2

FinalOutputFolder=$3

for fwdreads in $FilePathPrefix*val_1.fq.gz  ; do 

	fwdsuffix="_val_1.fq.gz"
	revsuffix="_val_2.fq.gz"
	revreads=${fwdreads/$fwdsuffix/$revsuffix}
	
	basestring=$(basename $fwdreads)
	idstring=${basestring/$fwdsuffix/""}
	
	outputfoldershovill=$Outputfolder$idstring

	shovillcontigs=$outputfoldershovill"/contigs.fa"
	

	newcontigplace=$FinalOutputFolder$idstring".fasta"


	shovill --outdir $outputfoldershovill --R1 $fwdreads --R2 $revreads
	cp $shovillcontigs $newcontigplace
done
