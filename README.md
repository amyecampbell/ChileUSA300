# ChileUSA300
Scripts used in USA300 SAE in Chile paper -- original versions of many of these scripts are in MRSA_Chile repo

## Assembly & QC

<br/>

bash [FastQC.sh](Assembly_QC/FastQC.sh) FastQCoutput rawreads

*Run with 279 raw read samples for isolates from Chile (in two batches, Nov '23 and Feb '24). Runs and outputs fastQC html reports.*

<br/>


bash [TrimReads.sh](Assembly_QC/TrimReads.sh) rawreads/ trimmedreads 1 2 

*Trim reads (10 bases off each end) using trim_galore, removing reads <50 bases in length*

<br/>

bash [Assemble.sh](Assembly_QC/Assemble.sh) trimmedreads ShovillOutput assemblies

*Assemble our 279 initial genomes from Chile (ST8 and non-ST8) using Shovill in paired read mode. Set up to give the file path prefix to all the samples you want to assemble (so you can break it up by prefix). Subsequently, copy these all into /scr1/users/campbela12/ST105/contigs_all_preQC/ as <genomeID>.fasta*

<br/>

bash [MashInfo.sh](Assembly_QC/MashInfo.sh)

*Make Mashinspected.tab, which gives some more taxonomic context to the hits. This will allow us to make sense of the non-aureus hits we get during MASH contamination check step-- e.g., is it mapping to S. haemolyticus because it's contaminated/the wrong species? Or is it mapping to a very closely related plasmid or phage in S. haemolyticus?*

<br/>

[RunMash.sh](Assembly_QC/RunMash.sh)
*Run mash screen with winner-take-all strategy with sketch database k=21, s=1000, RefSeq release 70. This yields output files for each assembly (/scr1/users/campbela12/ST105/MashOutputs_All/ or DataArchive/MashOutputs_All/)* 

<br/>

[Altogether_QC.R](Assembly_QC/Altogether_QC.R)

<br/>

*Put our 279 genomes from Chile (ST8 and non-ST8) and 62 StaphNET ST8 genomes through filters based on MASH contamination test, CheckM (<5% contamination, > 95% completeness), genome size (within 2 SDs of mean length in Genbank, 2863292). All StaphNET ST8-designated genomes pass, all but PP.3469, PP.3492, PP.3552, PP.3654 of the 279 Chilean MRSA genomes pass. PP.3646 was removed from this analysis due to mapping uncertainty.*

## Sequence typing

*Run MLST to assess % of each ST in 2022 isolates and to filter to ST8s and untyped SNVs of ST8*


## ST8 Phylogeny

[MakeSnippyInput.sh](TreeScripts/MakeSnippyInput.sh) SnippyInput.tab DataArchive/ST8Sequences/ 

<br/>

[MakeSnippyScript.sh](TreeScripts/MakeSnippyScript.sh)

<br/>

*Run Snippy on 433 ST8 genomes with GCA_000017085.1 as a reference genome.*

<br/>

[SnippyAlignment_ST8_Phylogeny.sh](TreeScripts/SnippyAlignment_ST8_Phylogeny.sh)

<br/>

Snippy "core" output files were moved to SnippyOutput/ for downstream steps

<br/>


*Run RaxML*
*Run ClonalFrameML*
*Run RaxML on recombination-masked phylogeny*
*Run Treeshrink (repeat previous steps after removing ERR9884500)*
*Run bootstrapping support*
*Run Allele typing script to type by Bianco et al. Alleles*

## SAE Tree: Phylogeny, Molecular Clock Analysis, Ancestral Reconstruction
*Run Snippy steps*
*Run RaxML*
*Run ClonalFrameML*
*Run RaxML on recombination-masked phylogeny*
*Run bootstrapping support*
*BactDating using tree*
*Ancestral recon with PastML*
*Data vis script for molecular clock and PastML results*

## COMER analysis 
*Run BLAST*
*Run BLAST*
*ReadMapAllGenomes.sh to align reads directly to SCCmec IVc --> COMER reference*
*R script to visualize read depth aligned to SCC mec through COMER*

## CURED
  *Run with 100/100 specificity/sensitivity cutoff*
 CURED_Main.py –case_control_file case_control.csv –extension fasta –genomes genomes/ 

 *Run with 0 specificity cutoff*
 CURED_Main.py –case_control_file case_control.csv –genomes genomes/ --sensitivity 100 –specificity 0 –extension fasta
