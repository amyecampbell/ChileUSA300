# ChileUSA300
Scripts used in USA300 SAE in Chile paper -- original versions of many of these scripts are in MRSA_Chile repo

## 0. Assembly, QC, basic stats

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
bash Assembly_QC/CheckMCall.sh /scr1/users/campbela12/ChileST108/Tree03_2024/Sequences/ /scr1/users/campbela12/ChileST108/Tree03_2024/CheckMoutput/ /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ChileST108/Tree03_2024/CheckM.txt

# StaphNET ST8 genomes
bash Assembly_QC/CheckMCall.sh /scr1/users/campbela12/ChileST108/StaphNET/contigs/ /scr1/users/campbela12/ChileST108/StaphNET/CheckM/ /scr1/users/campbela12/ST105/StaphAureusCheckMmarkers /scr1/users/campbela12/ChileST108/StaphNET/CheckMStaphNET_ST8.txt

```
*Run checkM with S aureus markers.*

<br/>

### MASH 

```
# RefSeqSketchesDefaults.msh is the k=21, s=1000 mash sketch database for RefSeq release 70 (255MB)
bash Assembly_QC/MashInfo.sh

grep '\.-p' /scr1/users/campbela12/downloadedDBs/Mashinspected.tab | sort -grk2 > MASHPlasmids_StartWithP.txt

grep 'unnamed-' /scr1/users/campbela12/downloadedDBs/Mashinspected.tab | sort -grk2 > MASHPlasmids_StartWithUnnamed.txt

awk '{ if ($2 < 200000) { print $3} }' MASHPlasmids_StartWithP.txt > MASHPlasmids_StartWithP_Less200kb.txt

awk '{ if ($2 < 200000) { print $3} }' MASHPlasmids_StartWithUnnamed.txt > MASHPlasmids_StartWithUnnamed_Less200kb.txt

cat MASHPlasmids_StartWithP_Less200kb.txt MASHPlasmids_StartWithUnnamed_Less200kb.txt >> MASH_likely_plasmids.txt

```


*Make Mashinspected.tab, which gives some more taxonomic context to the hits (tells us the full name of the sequence). This will allow us to make sense of the non-aureus hits we get during MASH contamination check step-- e.g., is it mapping to S. haemolyticus because it's contaminated/the wrong species? Or is it mapping to a very closely related plasmid in S. haemolyticus? Then, make a list of hit sources that are likely plasmids (start with the lowercase p and are less than 200kb, start with 'unnamed' and are less than 200kb (saved in DataArchive/MashFiles/MASH_likely_plasmids.txt)*

<br/>

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

<br/>

### TreeShrink

```
sh Assembly_QC//RunTreeShrink.sh

```

*Root a preliminary phylogeny and then use the resulting phylogeny to identify outlier branches. Note: this step was actually performed after the original run of the ST8 phylogeny scripts, and was used to identify that ERR9884473 had an outlier branch length(subsequently removed).*

<br/>

## 1. Figure 1 : Temporal analysis of MLST, antibiotic resistance, SAE over time

### Sequence typing

```
bash AlleleTyping/SequenceType.sh /scr1/users/campbela12/ST105/contigs_Feb/ /scr1/users/campbela12/ST105/ST105Tree_02_2024/ MLST_FebGenomes.txt

bash AlleleTyping/SequenceType.sh /scr1/users/campbela12/ST105/contigs/ /scr1/users/campbela12/ST105/

# Run later (after construction of ST8 phylogeny) to determine accuracy of CURED method
bash AlleleTyping/SequenceType.sh /scr1/users/campbela12/ChileST108/CURED/Contigs2023/SCLs/  /scr1/users/campbela12/ChileST108/CURED/ SCLmlst.txt 

```

*Run MLST to assess % of each ST in 2022 isolates and to filter to ST8s and untyped SNVs of ST8 (run separately on the three different sequencing batches). The following genomes were untyped by mlst, subsequently found via PubMLST searches (Feb 2024) to be SNVs and assigned to the following STs for downstream analysis purposes: PP.3354/SCL13027 -> ST5
PP.3515/SCL14812 -> ST8
PP.3555/SCL15367 -> ST105
PP.3589/SCL15881 -> ST105
PP.3598/SCL15958 -> ST5*

<br/>

### CURED to identify diagnostic restriction digests for USA300-SAE (via https://github.com/microbialARC/CURED)

```
CURED_Main.py –case_control_file case_control.csv –extension fasta –genomes genomes/ 
```

*Initially run CURED main script with default 100/100 specificity/sensitivity cutoff with DataArchive/CURED/local/case_control.csv input and did not find any 100% specific k-mers.*

```
 CURED_Main.py –case_control_file case_control.csv –genomes genomes/ --sensitivity 100 –specificity 0 –extension fasta
```
 *Rerun  with 0 specificity cutoff to identify non-specific k-mers. Using the k-mers of varying specificity listed in DataArchive/CURED/local/kmer_analysis/specificity_0/UniqueKmers.txt, filter report to identify those present in the fewest control genomes. Calculate corresponding specificity threshold (99%).*

```
CURED_Main.py –case_control_file case_control.csv –genomes genomes/ --sensitivity 100 –specificity 99 –extension fasta
```
*Rerun CURED with the updated specificity threshold(99%).*

<br/>

```
for file in $(cat DataArchive/CURED/local/kmer_analysis/specificity_0/UniqueKmers.txt); do grep -w "$file" DataArchive/CURED/local/kmer_analysis/specificity_99/temp_files/*.pyseer; done > DataArchive/CURED/local/specificity_99/grepped_output.txt

CURED/parse_grepped_op.py DataArchive/CURED/local/kmer_analysis/specificity_99/grepped_output.txt > DataArchive/CURED/local/kmer_analysis/specificity_99/parsed_grepped_output.txt

CURED/sort_cases_and_controls.py case_control.csv parsed_grepped_output.txt > DataArchive/CURED/local/kmer_analysis/specificity_99/controls_breaking_specificity.txt
```
<br/>


### Visualization
`
```
Rscript Fig1a-c_PctSTbyYear.R
```

*Produce figure 1A-C with % ST and resistance data by year (stored in DataArchive/Susceptibility_ST_Data_1999_2022.csv), SAE by year (based on ST8 tree-assigned clades stored in DataArchive/Figure2/TreeBasedCladesST8March_V2.csv and CURED estimate of 22.8% in 2023).*

<br/>

## 2. Figure 2: ST8 Phylogeny

### Snippy 4.6.0
```
# Same sequences are in DataArchive/ST8Sequences/
bash TreeScripts/MakeSnippyInput.sh Snippy_ST8_04_24.tab /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/Sequences/

bash TreeScripts/MakeSnippyScript.sh  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/GCA_000017085.1_reference.fasta Snippy_ST8_04_24.tab SnippyAlignment_ST8_Phylogeny.sh

# Manually add conda activation line, then run on cluster with 16G memory and 16 threads (sbatch -c 16 --mem=16G)
bash SnippyAlignment_ST8_Phylogeny.sh

# Moved output alignments to tree folder
mv *core*  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/

bash TreeScripts/SnippyClean.sh /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full.aln  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full_cleaned.aln

```

*Write a script that will run Snippy on 433 ST8 genomes with GCA_000017085.1 as a reference genome (SnippyAlignment_ST8_Phylogeny) and run it. Move the alignments output by snippy into a folder called ST8FullTree_04_2024, and then run the snippy-clean function to replace special characters with N for RaxML's use.*

<br/>

### RAxML 8.2.13 

```
# Provide job with 72 threads, 90G of memory via sbatch -c 72 --mem=90G
bash TreeScripts/RaxML_FirstRun.sh  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full_cleaned.aln 72 RaxML_ST8_April_PreCFML 100
```

*Run RaxML with 100 random starts and GTRGAMMA model*

<br/>

### ClonalFrameML 1.12 

```
# Provided with via 72G memory via sbatch --mem=72G
bash TreeScripts/RunCFML.sh /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/Trees/RAxML_bestTree.RaxML_ST8_April_PreCFML /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/core.full_cleaned.aln ST8AprilCFML

# Move all outputs of RaxML to the tree folder
mv *ST8AprilCFML* /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/
```
*Run clonalframeML on snippy alignment and ML phylogeny to mask recombinant regions. Outputs ST8AprilCFML.filtered.fasta*

<br/>

###  RAxML 8.2.13 following clonal frame ML 

```
# Provide with 72 threads and 90G memory via sbatch -c 72 ---mem=90G 
bash TreeScripts/RaxML_FirstRun.sh  /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/ST8AprilCFML.filtered.fasta 72 RaxML_ST8_April_PostCFML 100

```
*Run RaxML on recombination-masked phylogeny*

<br/>


```
# Provide with 64 threads and 96G of memory
TreeScripts/RaxML_Bootstrap.sh /scr1/users/campbela12/ChileST108/ST8FullTree_04_2024/ST8AprilCFML.filtered.fasta 64 RaxML_ST8_April_PostCFML RaxML_ST8_April_PostCFML_BootstrapTree RaxML_ST8_April_PostCFML_BootstrapParts 100

```
*Run bootstrapping support with 100 partitions*

<br/>

### Identify presence of [Bianco et al. (2023)'s](https://doi.org/10.3389/fcimb.2023.1081070) diagnostic alleles for USA300 clades in ST8 tree genomes

```
# 16 threads 16G memory
bash Bakta_ByFileList.sh SequenceList.txt /scr1/users/campbela12/ChileST108/Bakta/ /scr1/users/campbela12/ChileST108/GFFs/ /scr1/users/campbela12/ChileST108/ST8Genomes/
```
*Annotate genomes on ST8 tree*

<br/>


```
# ST8 genomes from 2022 collection
python3 Allele_Typing_BiancoEtAl.py

# Pre-2022 genomes from SCL collection
python3 Allele_Typing_BiancoEtAl_Pre2022.py

# ST8 genomes from StaphNET (10.1099/mgen.0.001020)
python3 Allele_Typing_BiancoEtAl_staphNET.py

# All genomes on ST8 tree (above 3 groups +)
python3 AlleleTyping/Allele_Typing_BiancoEtAl_AllST8s.py

```
*Run Allele typing script to type by Bianco et al. Alleles -- outputs presence (1) or absence(0) of each allele in each genome*

<br/>

```
Rscript Figure2_Table1_BiancoAlleleTesting.R
```
*Make mappings of leaf to Bianco et al. (2023) allele presence for ITOL visualization, calculate sensitivity/specificity of alleles.*
<br/>


## 3. Figure 3: SAE Phylogeny, molecular clock analysis, ancestral reconstruction 

<br/>

### Phylogeny

```
# Where Sequences/ contains the 154 sequences listed in "IsolatesUsed/Figures3_4_SAE Phylogeny" except for CA12.fasta, which is kept separate as the reference. 
bash TreeScripts/MakeSnippyInput.sh Snippy_ST8_SAE_V4.tab MolecularClockUTD_V2/Sequences/

bash TreeScripts/MakeSnippyScript.sh CA12.fasta Snippy_ST8_SAE_V4.tab MolecularClockTree/RunSnippy_ST8_SAE_V4.sh

# 16 threads and 20G memory
bash MolecularClockTree/RunSnippy_ST8_SAE_V4.sh

mv *core*  MolecularClockUTD_V2/

bash TreeScripts/SnippyClean.sh MolecularClockUTD_V2/core.full.aln  MolecularClockUTD_V2/core.full_cleaned.aln
```
*Run Snippy steps: make script based on list of 155 sequences, run script to generate alignment, clean alignment.*

<br/>

```
# 64 threads 64G memory
bash RaxML_FirstRun.sh MolecularClockUTD_V2/core.full_cleaned.aln 64 RaxML_SAE_April_V4 100
mv MolecularClockUTD_V2/*RaxML_SAE_April_V4* MolecularClockUTD_V2/Trees/
```
*Run RaxML before ClonalFrameML.*

<br/>


```

# 64G memory
bash TreeScripts/RunCFML.sh MolecularClockUTD_V2/Trees/RAxML_bestTree.RaxML_SAE_April_V4 MolecularClockUTD_V2/core.full_cleaned.aln SAE_CFML_V4
```
*Run ClonalFrameML to mask recombinant regions*

<br/>


```

# 64 threads and 64G memory
bash TreeScripts/RaxML_FirstRun.sh MolecularClockUTD_V2/SAE_CFML_V4.filtered.fasta 64 RaxML_SAE_April_V4_PostCFML 100

# 64 threads and 64G memory
bash TreeScripts/RaxML_Bootstrap.sh MolecularClockUTD_V2/SAE_CFML_V4.filtered.fasta 64 RaxML_SAE_April_V4_PostCFML RaxML_SAE_April_V4_PostCFML_BootstrapTree RaxML_SAE_April_V4_PostCFML_BootstrapParts 100
```
*Run RaxML on recombination-masked phylogeny, conduct bootstrapping*

<br/>

### Molecular clock analysis 

```
Rscript MolecularClockTree/BactDating.R
```
*BactDating using tree from RAxML. Outputs strictgammarun100million_RootFixed_April_V4.rda (R data object) and StrictClock100MillDatesFixed.newick (edited in FigTree v1.4.4 to produce strictClock100MillTree_StrippedDown.newick which replaces HPD with single branch length estimate)*

<br/>

### Ancestral reconstruction of geographic state
```
bash MolecularClockTree/PastML_Run_Full_MolecClock.sh 
```
*Ancestral recon with PastML using BactDating output tree. Outputs DataArchive/Figure 3/fullSAEpastMLoutputMolecClock/ which includes marginal probabilities file marginal_probabilities.character_Location.model_F81.tab needed for visualization*

<br/>

### Visualization of SAE subtree with branch lengths corresponding to time, internal nodes labeled by >50% MPP locations

```
Rscript Figure3_SAEsubtree.R
```
*Data vis script for molecular clock and PastML results*

<br/>

## Figure 4 -- Coverage of COMER and Copper Survival 


### COMER analysis 
*Run BLAST*
*Run BLAST*
*ReadMapAllGenomes.sh to align reads directly to SCCmec IVc --> COMER reference*
*R script to visualize read depth aligned to SCC mec through COMER*
