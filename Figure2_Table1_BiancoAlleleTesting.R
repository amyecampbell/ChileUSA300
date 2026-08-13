# Amy Campbell
# Figure 2: ITOL annotation tracks for % of different alleles present in full ST8 tree
# Table 1: Sensitivities and specificities of Bianco et al. (2023) alleles among Chilean ST8 and StaphNET ST8 isolates

library(dplyr) # version 1.1.4
library(stringr) # version 1.5.1



# Whole tree (done altogether for whole ST8 tree)
alleles_fulltree = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/AllelePresenceST8_FullST8Tree.csv")


# Just 2022 and pre-2022 Chilean, StaphNET as sensitivity/specificity test set
####################################################################################

# StaphNET ST8s
StaphNET_typing = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/AllelePresenceST8_StaphNET.csv")

# Pre-2022 Chilean ST8s
Pre2022Typing = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/AllelePresenceST8_Pre2022.csv")

# ST8s on 2022 isolates 
alleles2022 = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/AllelePresenceST8.csv")


#########################################
# Calculating sensitivity and specificity 
#########################################

cols_to_pull  = c("PPid","NAE_group_25850","NAE_YefM","NAE_PdhD","NAE_group_35132","NAE_group_1523","SAE_group_9725", "SAE_group_1496", 
                  "SAE_NAE_Der_2","SAE_NAE_DesR","SAE_NAE_NagE","SAE_NAE_SdcS","SAE_NAE_ComEC","USA300_group_3812", "USA300_AdhR","USA300_PchA",
                  "USA300_HisG","PEB1_HisB")

Pre2022Typed = Pre2022Typing[,cols_to_pull]
alleles2022Typed = alleles2022[,cols_to_pull]
StaphNET_typed = StaphNET_typing[,cols_to_pull]





Altogether_TestSet = rbind(StaphNET_typed, Pre2022Typed,alleles2022Typed)


TreeInfo = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/AnnotationCSVs/TreeClades.csv",header=T)
TreeInfo = read.csv("/Users/campbela12/Documents/PlanetGithubs/ChileUSA300/DataArchive/Figure\ 2/TreeBasedCladesST8March_V2.csv",header=F)
colnames(TreeInfo) = c("PPid", "clade")



Altogether_TestSet = Altogether_TestSet %>% left_join(TreeInfo,by="PPid")
Altogether_TestSet = Altogether_TestSet %>% filter(!(PPid %in% c("PP.3646","ERR9884473")))

NAE_cols = c("NAE_group_25850","NAE_YefM","NAE_PdhD","NAE_group_35132","NAE_group_1523")
SAE_cols = c("SAE_group_9725","SAE_group_1496")
SAENAE_cols = c("SAE_NAE_Der_2","SAE_NAE_DesR","SAE_NAE_NagE","SAE_NAE_SdcS","SAE_NAE_ComEC")
USA300_cols = c("USA300_group_3812","USA300_AdhR","USA300_PchA")
PEBcols = c("PEB1_HisB")

# Test NAE sensitivity specificity 
##################################
PositivesNAE = Altogether_TestSet %>% filter(clade=="NAE")
NegativesNAE = Altogether_TestSet %>% filter(clade!="NAE")

colSums(NegativesNAE[,NAE_cols])

NAEspecificity = 1-(colSums(NegativesNAE[,NAE_cols])/nrow(NegativesNAE))
NAEsensitivity = colSums(PositivesNAE[,NAE_cols])/nrow(PositivesNAE)
# NegativesNAE %>% filter(NAE_group_1523==1) PP.3341 contains the group_1523 (USA300HOU_1426) allele despite being an outgroup 

# Test SAE sensitivity specificity
##################################
PositivesSAE = Altogether_TestSet %>% filter(clade=="SAE")
NegativesSAE = Altogether_TestSet %>% filter(clade!="SAE")

SAEspecificity = 1-(colSums(NegativesSAE[,SAE_cols])/nrow(NegativesSAE))
SAEsensitivity = colSums(PositivesSAE[,SAE_cols])/nrow(PositivesSAE)

# SAE/NAE combo sensitity and specificity 
#########################################
PositivesSAENAE = Altogether_TestSet %>% filter(clade %in% c("NAE", "SAE"))
NegativesSAENAE = Altogether_TestSet %>% filter(!(clade %in% c("NAE", "SAE")))

SAENAEspecificity = 1-(colSums(NegativesSAENAE[,SAENAE_cols])/nrow(NegativesSAENAE))
SAENAEsensitivity = colSums(PositivesSAENAE[,SAENAE_cols])/nrow(PositivesSAENAE)

# PEB1 sensitivity and specificity 
##################################
PEBPositives = Altogether_TestSet %>% filter(clade=="PEB1")
PEBNegatives = Altogether_TestSet %>% filter(clade!="PEB1")
PEB1specificity = 1-(sum(PEBNegatives[,"PEB1_HisB"])/nrow(PEBNegatives))

PEB1sensitivity = sum(PEBPositives[,"PEB1_HisB"])/nrow(PEBPositives)

PEB1sensitivity = c(PEB1sensitivity)
names(PEB1sensitivity) = c("PEB1_HisB")


PEB1specificity = c(PEB1specificity)
names(PEB1specificity) = c("PEB1_HisB")


# USA300 sensitivity and specificity 
####################################

PositivesUSA300 = Altogether_TestSet %>% filter(clade!="Outgroup")
NegativesUSA300 = Altogether_TestSet %>% filter(clade=="Outgroup")

USA300specificity = 1-(colSums(NegativesUSA300[,USA300_cols])/nrow(NegativesUSA300))
USA300sensitivity = colSums(PositivesUSA300[,USA300_cols])/nrow(PositivesUSA300)

# Combine
##########
sensitivities = c(NAEsensitivity, SAEsensitivity, SAENAEsensitivity,PEB1sensitivity, USA300sensitivity)
specificities = c(NAEspecificity, SAEspecificity, SAENAEspecificity, PEB1specificity, USA300specificity)

sens_spec_DF = data.frame(Sensitivity = sensitivities, Specificity = specificities)
sens_spec_DF$Specificity = 100*sens_spec_DF$Specificity
sens_spec_DF$Sensitivity = 100*sens_spec_DF$Sensitivity


write.csv(sens_spec_DF, "FigureOutput/Table1.csv")


#####################################################################################
# Getting % each Bianco et al. allele present for ITOL visualization of full ST8 tree
#####################################################################################
alleles_fulltree$PctNAE = 100*(rowSums(alleles_fulltree[c("NAE_group_25850", "NAE_YefM", "NAE_PdhD","NAE_group_35132", "NAE_group_1523")])/5)
alleles_fulltree$PctSAE = 100*(rowSums(alleles_fulltree[c("SAE_group_9725", "SAE_group_1496")])/2)
alleles_fulltree$PctNAE_SAE = 100*(rowSums(alleles_fulltree[c("SAE_NAE_Der_2", "SAE_NAE_DesR", "SAE_NAE_NagE","SAE_NAE_SdcS","SAE_NAE_ComEC")])/5)
alleles_fulltree$PctUSA300 = 100*(rowSums(alleles_fulltree[c("USA300_group_3812", "USA300_group_3812", "USA300_PchA","USA300_HisG")])/4)
alleles_fulltree$PctPEB1 = alleles_fulltree$PEB1_HisB*100
Pcts = (alleles_fulltree %>% select(PPid, PctUSA300, PctPEB1, PctNAE, PctSAE ))

write.table(Pcts, "~/Documents/Planet/ST105/ST8/AnnotationCSVs/BiancoEtAlAlleles.csv", sep=",",quote=F,row.names=F,col.names=F)

write.table(Pcts %>% select(PPid, PctUSA300), "~/Documents/Planet/ST105/ST8/AnnotationCSVs/BiancoEtAllUSA300.csv", sep=",",quote=F,row.names=F,col.names=F)
write.table(Pcts %>% select(PPid, PctPEB1), "~/Documents/Planet/ST105/ST8/AnnotationCSVs/BiancoEtAllPeb1.csv", sep=",",quote=F,row.names=F,col.names=F)
write.table(Pcts %>% select(PPid, PctNAE), "~/Documents/Planet/ST105/ST8/AnnotationCSVs/BiancoEtAllNAE.csv", sep=",",quote=F,row.names=F,col.names=F)
write.table(Pcts %>% select(PPid, PctSAE), "~/Documents/Planet/ST105/ST8/AnnotationCSVs/BiancoEtAllSAE.csv", sep=",",quote=F,row.names=F,col.names=F)


TreeBasedAssignments = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/Trees/TreeBasedCladesST8March_V2.csv",header=F)
TreeBasedAssignments$V1 = sapply(TreeBasedAssignments$V1, function(x) str_replace_all(x, " ","_"))
TreeBasedAssignments$PPid = TreeBasedAssignments$V1

alleles_fulltree = alleles_fulltree %>% left_join(TreeBasedAssignments, by="PPid")

write.csv(alleles_fulltree,"FigureOutput/AnnotationMetadata.csv")
