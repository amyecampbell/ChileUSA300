#!/bin/bash

#conda activate blastenv

outputprefix="./ST8/BlastST8s/"
mkdir -p $outputprefix

#filename=$(basename $testfile)
#ppID=${filename/".fasta"/""}

blastdbpath="./ST8/BlastDBs/"

comerpath="DataArchive/Figure\ 4/CA12_COMER.fasta"
mkdir -p $blastdbpath

export BLASTDB=$blastdbpath

for testfile in ./DataArchive/Sequences/*.fasta; do

	filename=$(basename $testfile)
	ppID=${filename/".fasta"/""}
	echo $ppID
	makeblastdb -in $testfile -dbtype nucl -out $ppID

	mv $ppID* $blastdbpath/

	blastn -db $ppID -query $comerpath -out $outputprefix$ppID"_comer.tab" -outfmt "6 qseqid slen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
	myprefix=$ppID"\t"
	#sed 's/.*/$ppID"\t"/' $outputprefix$ppID"_comer.tab"
	#perl -pe 's/^/$ppID\t/' $outputprefix$ppID"_comer.tab" 
	awk -v prefix=$myprefix '{print prefix $0}' $outputprefix$ppID"_comer.tab" > $outputprefix$ppID"_comer_mod.tab"
	#sed 's/^/$prefix_/' $outputprefix$ppID"_comer.tab" > $outputprefix$ppID"_comer_mod.tab" 
	#awk '{print $prefix $0}' $outputprefix$ppID"_comer.tab" 

done
cat $outputprefix*comer_mod.tab > $outputprefix"COMERall.tab"

mv  $outputprefix"COMERall.tab" DataArchive/Figure\ 4/
