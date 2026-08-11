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

*Make Mashinspected.tab, which gives some more taxonomic context to the hits (tells us the full name of the sequence). This will allow us to make sense of the non-aureus hits we get during MASH contamination check step-- e.g., is it mapping to S. haemolyticus because it's contaminated/the wrong species? Or is it mapping to a very closely related plasmid in S. haemolyticus? Then, make a list of hit sources that are likely plasmids (start with the lowercase p and are less than 200kb, start with 'unnamed' and are less than 200kb (saved in DataArchive/MashFiles/MASH_likely_plasmids.txt)*

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

### TreeShrink

```
sh Assembly_QC//RunTreeShrink.sh

```

*Root a preliminary phylogeny and then use the resulting phylogeny to identify outlier branches. Note: this step was actually performed after the original run of the ST8 phylogeny scripts, and was used to identify that ERR9884473 had an outlier branch length(subsequently removed).*

### Sequence typing


```
bash TreeScripts/SequenceType.sh /scr1/users/campbela12/ST105/contigs_Feb/ /scr1/users/campbela12/ST105/ST105Tree_02_2024/ MLST_FebGenomes.txt

bash TreeScripts/SequenceType.sh /scr1/users/campbela12/ST105/contigs/ /scr1/users/campbela12/ST105/

# Run later (after construction of ST8 phylogeny) to determine accuracy of CURED method
bash TreeScripts/SequenceType.sh /scr1/users/campbela12/ChileST108/CURED/Contigs2023/SCLs/  /scr1/users/campbela12/ChileST108/CURED/ SCLmlst.txt 

```

*Run MLST to assess % of each ST in 2022 isolates and to filter to ST8s and untyped SNVs of ST8 (run separately on the three different sequencing batches). The following genomes were untyped by mlst, subsequently found via PubMLST searches (Feb 2024) to be SNVs and assigned to the following STs for downstream analysis purposes: PP.3354/SCL13027 -> ST5
PP.3515/SCL14812 -> ST8
PP.3555/SCL15367 -> ST105
PP.3589/SCL15881 -> ST105
PP.3598/SCL15958 -> ST5*


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

### Visualization

```
Rscript Fig1a-c_PctSTbyYear.R
```

*Produce figure 1A-C with % ST and resistance data by year (stored in DataArchive/Susceptibility_ST_Data_1999_2022.csv), SAE by year (based on ST8 tree-assigned clades stored in DataArchive/Figure2/TreeBasedCladesST8March_V2.csv and CURED estimate of 22.8% in 2023).*

## Figure 2 -- ST8 Phylogeny

```

# Same sequences are in DataArchive/ST8Sequences/
bash TreeScripts/MakeSnippyInput.sh Snippy_ST8_04_24.tab /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/Sequences/

bash TreeScripts/MakeSnippyScript.sh  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/GCA_000017085.1_reference.fasta Snippy_ST8_04_24.tab SnippyAlignment_ST8_Phylogeny.sh

# Run on cluster with 16G memory and 16 threads (sbatch -c 16 --mem=16G)
bash SnippyAlignment_ST8_Phylogeny.sh

# Moved output alignments to tree folder
mv *core*  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/

sbatch /home/campbela12/Documents/MRSA_Chile/Trees/SnippyClean.sh /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full.aln  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full_cleaned.aln

```

<br/>

*Write a script that will run Snippy on 433 ST8 genomes with GCA_000017085.1 as a reference genome. Outputs SnippyAlignment_ST8_Phylogeny.sh (originally run with path /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/Sequences/ instead of DataArchive/ST8Sequences/) *


*Run Snippy alignment with GCA_000017085.1 as a reference genome. Following this Snippy "core" output files were moved to SnippyOutput/ for downstream steps*

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
