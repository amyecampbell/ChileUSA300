# Amy Campbell
# 2024
# Quality control check on all 279 sequenced chilean MRSA genomes together
library(stringr) # 1.5.1 
library(dplyr) # 1.1.4

# Quast (version 5.2.0)
##########################
Quast = read.csv2("Documents/Planet/ST105/QC_files_all279/Quast/transposed_report.tsv",sep="\t",header=T)

# within [(mean S. aureus genome size) +/- 2x(SD S. aureus genome size)] based on Sean's Feb 20, 2024 query of Genbank with metadata
Quast$PassedSize = if_else(Quast$Total.length >= (2863292 - 2*121772) & Quast$Total.length <= (2863292 + 2*121772), "Yes","No")

# CheckM (version 1.2.2)
##########################
checkM_Nov = read.table("Documents/Planet/ST105/QC_files_all279/checkMoutput_Nov23.txt",skip=10,header=F)
checkM_Feb = read.table("Documents/Planet/ST105/QC_files_all279/checkMoutput_Feb24.txt",skip=10,header=F)
CheckM_all = (rbind(checkM_Nov, checkM_Feb))
CheckM_all$Assembly = CheckM_all$V1

CheckM_all$PassedCheckM = if_else((CheckM_all$V14 >= 95) & (CheckM_all$V15< 5.0), "Yes", "No" )

QC_Passes = Quast %>% select(Assembly, PassedSize) %>% left_join(CheckM_all %>% select(Assembly, PassedCheckM),by="Assembly")

# these checks are congruent
QC_Passes %>% filter(PassedSize==PassedCheckM)

# Mash (version 2.3) (looking for non-plasmids with >100/1000 hits and >95% identity)
#######################################################################################
MashResults = read.csv2("Documents/Planet/ST105/QC_files_all279/MASH/MRSAChileAll_ExcludeLikelPlasmids_annotated.tsv",sep="\t")
MashResults$Assembly = sapply(MashResults$File, function(x) str_remove_all(x, "_Mash.tab"))
MashResults$PassedMash= MashResults$Passed
QC_Passes = QC_Passes %>% left_join(MashResults %>% select(Assembly, PassedMash), by="Assembly")

Failed = QC_Passes %>% filter(PassedSize=="No" |  PassedCheckM=="No" | PassedMash=="No")
Passed = setdiff(QC_Passes$Assembly, Failed$Assembly)
write.table(Passed, "Documents/Planet/ST105/DataSummaries/PassedQC.txt",quote=F, col.names=F, row.names=F)

# Failed QC : PP.3469, PP.3492, PP.3552, PP.3654


# Check these metrics for the staphnet genomes
###############################################
CheckM_StaphNET = read.table("/Users/campbela12/Documents/Planet/ST105/ST8/ST8_StaphNET-SA-First-Survey/CheckMStaphNET_ST8.txt",skip=10,header=F)
MashStaphNET = read.csv2("/Users/campbela12/Documents/Planet/ST105/ST8/ST8_StaphNET-SA-First-Survey/ST8_staphNETMashResults.txt",sep="\t")
Quast_StaphNET = read.csv2("/Users/campbela12/Documents/Planet/ST105/ST8/ST8_StaphNET-SA-First-Survey/Quast_StaphNET.tsv",sep="\t",header=T)
Quast_StaphNET$PassedSize = if_else(Quast_StaphNET$Total.length >= (2863292 - 2*121772) & Quast_StaphNET$Total.length <= (2863292 + 2*121772), "Yes","No")

CheckM_StaphNET$PassedCheckM = if_else((CheckM_StaphNET$V14 >= 95) & (CheckM_StaphNET$V15< 5.0), "Yes", "No" )
table(MashStaphNET$PassedMASH)
table(Quast_StaphNET$PassedSize)
table(CheckM_StaphNET$PassedCheckM)

# all StaphNET ST8 genomes pass all 3 quality criteria 
