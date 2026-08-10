# ChileUSA300
Scripts used in USA300 SAE in Chile paper -- original versions of many of these scripts are in MRSA_Chile repo

## Assembly, QC, basic stats

<br/>

```
bash Assembly_QC/FastQC.sh FastQCoutput rawreads
```

*Run FastQC with 279 raw read samples for isolates from Chile (in two batches, Nov '23 and Feb '24). Runs and outputs fastQC html reports.*

<br/>

```
bash Assembly_QC/TrimReads.sh rawreads/ trimmedreads 1 2 
```

*Trim reads (10 bases off each end) using trim_galore, removing reads <50 bases in length*

<br/>

```
bash Assembly_QC/Assemble.sh trimmedreads ShovillOutput assemblies
```

*Assemble our 279 initial genomes from Chile (ST8 and non-ST8) using Shovill in paired read mode. Set up to give the file path prefix to all the samples you want to assemble (so you can break it up by prefix). Subsequently, copied these all into /scr1/users/campbela12/ST105/contigs_all_preQC/ as <genomeID>.fasta*

<br/>

### CheckM with S. aureus markers set (<5% contamination, >95% completeness)

```

# First run of 2022 isolates
bash Assembly_QC/CheckMCall.sh /scr1/users/campbela12/ST105/contigs/ /scr1/users/campbela12/ST105/CheckMoutput /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ST105/checkMoutput_Nov.txt

# Second run of 2022 isolates
bash Assembly_QC/CheckMCall.sh /scr1/users/campbela12/ST105/contigs_Feb/ /scr1/users/campbela12/ST105/CheckMoutputFeb/ /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ST105/CheckM_Feb24

# ST8 genomes (non-StaphNET)
bash CheckMCall.sh /scr1/users/campbela12/ChileST108/Tree03_2024/Sequences/ /scr1/users/campbela12/ChileST108/Tree03_2024/CheckMoutput/ /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ChileST108/Tree03_2024/CheckM.txt

# StaphNET ST8 genomes
bash /home/campbela12/Documents/MRSA_Chile/assembly/CheckMCall.sh /scr1/users/campbela12/ChileST108/StaphNET/contigs/ /scr1/users/campbela12/ChileST108/StaphNET/CheckM/ /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ChileST108/StaphNET/CheckMStaphNET_ST8.txt


```
### MASH 

```
bash Assembly_QC/MashInfo.sh
grep '\.-p' /scr1/users/campbela12/downloadedDBs/Mashinspected.tab | sort -grk2 > MASHPlasmids_StartWithP.txt

grep 'unnamed-' /scr1/users/campbela12/downloadedDBs/Mashinspected.tab | sort -grk2 > MASHPlasmids_StartWithUnnamed.txt

awk '{ if ($2 < 200000) { print $3} }' MASHPlasmids_StartWithP.txt > MASHPlasmids_StartWithP_Less200kb.txt

awk '{ if ($2 < 200000) { print $3} }' MASHPlasmids_StartWithUnnamed.txt > MASHPlasmids_StartWithUnnamed_Less200kb.txt

cat MASHPlasmids_StartWithP_Less200kb.txt MASHPlasmids_StartWithUnnamed_Less200kb.txt >> MASH_likely_plasmids.txt

<br/>

```

*Make Mashinspected.tab, which gives some more taxonomic context to the hits. This will allow us to make sense of the non-aureus hits we get during MASH contamination check step-- e.g., is it mapping to S. haemolyticus because it's contaminated/the wrong species? Or is it mapping to a very closely related plasmid in S. haemolyticus? Then, make a list of hit sources that are likely plasmids (start with the lowercase p and are less than 200kb, start with 'unnamed' and are less than 200kb (saved in DataArchive/MashFiles/MASH_likely_plasmids.txt)*

```
bash Assembly_QC/RunMash.sh
```

*Run mash screen with winner-take-all strategy with sketch database k=21, s=1000, RefSeq release 70. This yields output files for each assembly (/scr1/users/campbela12/ST105/MashOutputs_All/, stored in DataArchive/MashOutputs_All/)* 


<br/>

```
python3 Assembly_QC/Mash_contamination_checker_ModifiedAEC.py -i .95 -s 100 -p ST8_Genomes_Tree "Staphylococcus_aureus" --exclude_names_file /Users/campbela12/Documents/Planet/ST105/QC_files_all279/likelyPlasmids_Exclude_MASH.txt  /Users/campbela12/Documents/Planet/ST105/ST8/ST8_QC/MASHfullTree/

python3 Assembly_QC/Mash_contamination_checker_ModifiedAEC.py -i .95 -s 100 -p MRSAChileAll_ExcludeLikelPlasmids "Staphylococcus_aureus" --exclude_names_file /Users/campbela12/Documents/Planet/ST105/QC_files_all279/MASH/MASH_likely_plasmids.txt /Users/campbela12/Documents/Planet/ST105/QC_files_all279/MASH/MashOutputs_All

python3 Assembly_QC/Mash_contamination_checker_ModifiedAEC.py -i .95 -s 100 -p ST8_staphNET "Staphylococcus_aureus" --exclude_names_file /Users/campbela12/Documents/Planet/ST105/QC_files_all279/MASH/MASH_likely_plasmids.txt  /Users/campbela12/Documents/Planet/ST105/ST8/ST8_StaphNET-SA-First-Survey/MASHStaphNET
```

*Following mash screen, run contamination checker (A. Moustafa's script, modified by AEC), excluding likely plasmids on the 279 2022 genomes, the StapNET genomes, and all genomes included in the ST8 tree. After check for MASH contamination, CheckM contamination and completeness (<5% contamination, > 95% completeness), genome size (within 2 SDs of mean length in Genbank, 2863292), the following genomes which had initially been included were removed from analysis:
ERR134761 (Public ST8),
ERR134840 (Public ST8),
ERR137917 (Public ST8),
ERR137937 (Public ST8),
SRR497495 (Public ST8),
SCL14023/PP.3469 (2022 genome),
SCL14160/PP.3492 (2022 genome),
SCL15359/PP.3552 (2022 genome),
SCL16696/PP.3654  (2022 genome),
SCL16471/PP.3646 (removed for metadata mapping uncertainty)* 

### TreeNET

*Note: this step was actually performed after the original run of the ST8 phylogeny scripts, and was used to identify *

### Sequence typing

```

```

*Run MLST to assess % of each ST in 2022 isolates and to filter to ST8s and untyped SNVs of ST8*


## Figure 1 : Temporal analysis of MLST, antibiotic resistance, SAE over time

### CURED to identify diagnostic restriction digests for USA300-SAE (via https://github.com/microbialARC/CURED)

```
 CURED_Main.py –case_control_file case_control.csv –extension fasta –genomes genomes/ 
```

* Initially ran CURED main script with default 100/100 specificity/sensitivity cutoff, and did not find any 100% specific k-mers.*

```
 CURED_Main.py –case_control_file case_control.csv –genomes genomes/ --sensitivity 100 –specificity 0 –extension fasta
```

 *Re-ran with 0 specificity cutoff to identify non-specific k-mers*
```
CURED_Main.py –case_control_file case_control.csv –genomes genomes/ --sensitivity 100 –specificity 0 –extension fasta
```
*Additional CURED steps here*

### 



## Figure 2 -- ST8 Phylogeny

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


## Figure 3 -- SAE Phylogeny, molecular clock analysis, ancestral reconstruction 

### SAE Tree: Phylogeny, Molecular Clock Analysis, Ancestral Reconstruction

*Run Snippy steps*
*Run RaxML*
*Run ClonalFrameML*
*Run RaxML on recombination-masked phylogeny*
*Run bootstrapping support*
*BactDating using tree*
*Ancestral recon with PastML*
*Data vis script for molecular clock and PastML results*

## Figure 4 -- Coverage of COMER and Copper Survival 

### COMER analysis 
*Run BLAST*
*Run BLAST*
*ReadMapAllGenomes.sh to align reads directly to SCCmec IVc --> COMER reference*
*R script to visualize read depth aligned to SCC mec through COMER*
