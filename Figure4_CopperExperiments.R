# Amy Campbell
# 2026 (Some code is extracted from 2025 code in MRSA_Chile repo)
# Plotting 
library(tidyverse)
library(ggpubr)
library(stringr)
library(reshape2)


###################
# COPPER IN WATER 
##################

# Formatted the traditional way (one row per individual observation)
CopperData <- read.csv("DataArchive/RawData_CopperInWater.csv")

# Averaging two reps of the experiment I did back in 2024 since they were from the same colony
AmyExperimentCollapsed = CopperData %>% filter(Experiment == "Amy") %>% group_by(Copper, Strain) %>% summarize(Experiment="Amy",CFUml = mean(CFUml))

CopperDataFiltered = CopperData %>% filter(Experiment != "Amy")
CopperDataFiltered$log10cfu=log10(CopperDataFiltered$CFUml+1)

CopperData$log10cfu = log10(CopperData$CFUml + 1)

CopperDataFilteredMeans = CopperDataFiltered %>% group_by(Copper, Strain) %>% summarize(CFUml= mean(CFUml), meanLogCFUml= mean(log10cfu), sdlog10 = sd(log10cfu), n=n())

CopperDataFilteredMeans$selog10 = CopperDataFilteredMeans$sdlog10/sqrt(CopperDataFilteredMeans$n)

CopperDataRenamed = CopperDataFiltered
CopperDataRenamed = CopperDataRenamed %>% mutate(Strain = case_when(Strain=="M121" ~ "M121",
                                                                    Strain=="KO" ~"M121-delta-cop",
                                                                    Strain=="Comp3970" ~ "M121-delta-cop::cop",
                                                                    Strain=="EmptyVector" ~ "Vector"))
CopperDataFilteredMeansRenamed = CopperDataFilteredMeans %>% mutate(Strain = case_when(Strain=="M121" ~ "M121",
                                                                                       Strain=="KO" ~"M121-delta-cop",
                                                                                       Strain=="Comp3970" ~ "M121-delta-cop::cop",
                                                                                       Strain=="EmptyVector" ~ "Vector"))

CopperDataFilteredMeansRenamed$CFUml <- NULL

CopperDataRenamed <- CopperDataRenamed %>% left_join(CopperDataFilteredMeansRenamed, by=c("Copper","Strain"))

# Boxplot by strain
#####################
boxplot_copper = ggplot(CopperDataRenamed, aes(x=factor(Copper), y=log10cfu, fill=factor(Strain))) + geom_boxplot(alpha=.7) +  geom_jitter(data= CopperDataRenamed,aes(x=factor(Copper),y=log10cfu),
                                                                                                                                           position=position_jitterdodge(dodge.width =.75, jitter.width=.3)) + scale_fill_brewer(palette="Dark2")  + theme_classic()



ggsave(boxplot_copper + theme_classic(), file="FigureOutput/Fig4_BoxplotCopperInWater.pdf",width=8, height=5)



# Pairing is only possible for mutant vs. WT because early experiments only had mutant and WT
CopperDataRenamed10 = CopperDataRenamed %>% filter(Copper=="10µM")

attach(CopperDataRenamed10)
pairwise.wilcox.test(log10cfu, Strain, p.adjust.method = "BH")
detach(CopperDataRenamed10)

CopperDataRenamed0 = CopperDataRenamed %>% filter(Copper=="0µM")
attach(CopperDataRenamed0)
pairwise.wilcox.test(log10cfu, Strain, p.adjust.method = "BH")
detach(CopperDataRenamed0)



################
# COPPER IN RPMI 
################

Cfus_ml_24 = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/RPMIExperiments_2026/RPMI_CFUmls_24.csv")
Cfus_ml_0 = read.csv("/Users/campbela12/Documents/Planet/ST105/ST8/RPMIExperiments_2026/RPMIInitialCFUs.csv")

Cfus_ml_24$Strain = sapply(Cfus_ml_24$Strain.Concentration, function(x) str_split(x, "\\-")[[1]][1])

Cfus_ml_24$Concentration = sapply(Cfus_ml_24$Strain.Concentration, function(x) as.numeric(as.character(str_split(x, "\\-")[[1]][2])))

Renamed = Cfus_ml_24 %>% mutate(StrainName = case_when(Strain=="821" ~ "M121",
                                                       Strain=="1066" ~"M121-delta-cop",
                                                       Strain=="3970" ~ "M121-delta-cop::cop",
                                                       Strain=="3976" ~ "Vector"))


Cfus_ml_24 = Renamed
Cfus_ml_24$Strain_Exp = paste0(Cfus_ml_24$Strain, "_", Cfus_ml_24$Experiment)
Cfus_ml_0$Strain_Exp = paste0(Cfus_ml_0$Strain, "_", Cfus_ml_0$Experiment)

Cfus_ml_0_Merge = Cfus_ml_0 %>% select(Strain_Exp, Log10.CFU.mL)

colnames(Cfus_ml_0_Merge) = c("Strain_Exp", "log10_0h_pseudocount")

Cfus_ml_24 = Cfus_ml_24 %>% left_join(Cfus_ml_0_Merge, by="Strain_Exp")

Cfus_ml_24$LogFoldChange = Cfus_ml_24$log10plus1_CFUml_24H - Cfus_ml_24$log10_0h_pseudocount


Cfus_ml_24_Means = Cfus_ml_24 %>% group_by(Strain, Concentration) %>% summarize(StrainName=StrainName, Concentration=Concentration, avgLogChange = mean(LogFoldChange), sdfc = sd(LogFoldChange))


Cfus_ml_24Pt5 = Cfus_ml_24 %>% filter(Concentration %in% c(0, .5))
Cfus_ml_24_MeansPt5 = Cfus_ml_24_Means %>% filter(Concentration %in% c(0, .5))

FoldChangePlot0pt5 <- ggplot(data = Cfus_ml_24_MeansPt5, aes(x = factor(Concentration), y = avgLogChange, fill = factor(StrainName)) ,alpha=.7) +
  geom_bar(stat = "identity", position = position_dodge(.9)) + scale_y_continuous(breaks=seq(-8 , 2, 1)) +
  geom_errorbar(aes(ymin=avgLogChange-sdfc, ymax=avgLogChange+sdfc), position=position_dodge(.9), width=.25) + scale_fill_brewer(palette = "Dark2") + theme_classic() +
  geom_jitter(data= Cfus_ml_24Pt5,aes(x=factor(Concentration),y=LogFoldChange),
              position=position_jitterdodge(dodge.width =.75, jitter.width=.2))
  
 ggsave(FoldChangePlot0pt5, file="FigureOutput/Figure4_RPMIGrowth.pdf",width=7,height=5)
  
 
Cfus_ml_24_0 = Cfus_ml_24 %>% filter(Concentration == 0)
shapiro.test(Cfus_ml_24_0$LogFoldChange)
attach(Cfus_ml_24_0)
pairwise.t.test(LogFoldChange, Strain, p.adjust.method = "BH")
detach(Cfus_ml_24_0)


Cfus_ml_24_5 = Cfus_ml_24 %>% filter(Concentration == .5)
shapiro.test(Cfus_ml_24_5$LogFoldChange)
attach(Cfus_ml_24_5)
pairwise.t.test(LogFoldChange, Strain, p.adjust.method = "BH")
detach(Cfus_ml_24_5)


# Supplements
##############
FoldChangePlot <- ggplot(data = Cfus_ml_24_Means, aes(x = factor(Concentration), y = avgLogChange, fill = factor(StrainName)) ,alpha=.7) +
  geom_bar(stat = "identity", position = position_dodge(.9)) + scale_y_continuous(breaks=seq(-8 , 2, 1)) +
  geom_errorbar(aes(ymin=avgLogChange-sdfc, ymax=avgLogChange+sdfc), position=position_dodge(.9), width=.25) + scale_fill_brewer(palette = "Dark2") + theme_classic() +
  geom_jitter(data= Cfus_ml_24,aes(x=factor(Concentration),y=LogFoldChange),
              position=position_jitterdodge(dodge.width =.75, jitter.width=.2))

ggsave(FoldChangePlot, file="FigureOutput/Figure4supplement_RPMIGrowth_allconcentrations.pdf",width=12,height=8)

