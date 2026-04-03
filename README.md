# ChileUSA300
Scripts used in USA300 SAE in Chile paper -- original versions of many of these scripts are in MRSA_Chile repo

## Assembly & QC

<br/>

[FastQC.sh](Assembly_QC/FastQC.sh)

*Run with 279 raw read samples for isolates from Chile (in two batches, Nov '23 and Feb '24). Runs and outputs fastQC html reports.*

<br/>


[TrimReads.sh](Assembly_QC/TrimReads.sh)

*Trim reads (10 bases off each end) using trim_galore, removing *

<br/>

[Assemble.sh](Assembly_QC/Assemble.sh)

*Assemble our 279 initial genomes from Chile (ST8 and non-ST8) using Shovill in paired read mode. Set up to give the file path prefix to all the samples you want to assemble (so you can break it up by prefix)*

<br/>



[Altogether_QC.R](Assembly_QC/Altogether_QC.R)

*Put our 279 genomes from Chile (ST8 and non-ST8) and 62 StaphNET ST8 genomes through filters based on MASH contamination test, CheckM (<5% contamination, > 95% completeness), genome size (within 2 SDs of mean length in Genbank, 2863292). All StaphNET ST8-designated genomes pass, all but PP.3469, PP.3492, PP.3552, PP.3654 of the 279 Chilean MRSA genomes pass.*
