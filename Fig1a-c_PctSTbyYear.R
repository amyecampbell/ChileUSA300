# Amy Campbell
# Making plot of ST % by year in Chile invasive infections
library(dplyr) # version 1.1.4 
library(stringr) # version 1.5.1
library(RColorBrewer) # version 1.1-3
library(ggplot2) # version 3.4.4 
library(reshape2) # version 1.4.4
library(tidyr) # version 1.3.1
library(collapse) # version 2.0.12
library(fmsb) # version 0.7.6
library(ggstream) # version 0.1.0
library(ztable) # version0.2.3

# Read in data
##############
green="#4DAF4A"
red="#E31A1C"
blue="#386CB0"
pink="#E7298A"
turquoise="#66C2A5"
orange="#FF7F00"
brown="#A65628"
gold="#E6AB02"
purp="#984EA3"
yellow="#FFFF33"
lightpurp= "#BEAED4"
lightyellow="#FFF2AE"
gray = "#999999"
lightgreen="#B3DE69"
lightblue="#A6CEE3"
lightorange="#FDB462"

AllIsolates = read.csv("./DataArchive/Figure 1/Susceptibility_ST_Data_1999_2022.csv")
AllIsolatesDedup = AllIsolates %>% filter(IncludeDeduplicated == 1)


# Convert single nucleotide variants to their closest, and make an 'Other STs' category just for visualization
AllIsolates$ST = if_else(AllIsolates$ST=="~105", "105",AllIsolates$ST)
AllIsolates$ST = if_else(AllIsolates$ST=="~1472", "1472",AllIsolates$ST)
AllIsolates$ST = if_else(AllIsolates$ST %in% c("no descrito","-" ), "Untyped",AllIsolates$ST)

OtherSTs = names(table(AllIsolatesDedup$ST)[table(AllIsolatesDedup$ST) <=2])
OtherSTs = setdiff(OtherSTs, "Untyped")
AllIsolates$ST = if_else(AllIsolates$ST %in% OtherSTs, "Other STs",AllIsolates$ST)



AllIsolates = AllIsolates %>% mutate(colorscheme=case_when(ST=="5" ~blue,
                                                                     ST=="105" ~green ,
                                                                     ST=="8" ~ orange,
                                                                     ST=="72" ~yellow,
                                                                     ST=="1472" ~ turquoise,
                                                                     ST=="225"~ gold,
                                                                     ST== "30" ~ purp,
                                                                     ST == "239" ~ lightyellow,
                                                                     ST=="923"~ brown,
                                                                     ST=="2802"~ pink,
                                                                     ST=="97"~ red,
                                                                     ST == "Other STs"  ~ "#454545",
                                                                     ST=="Untyped" ~ "#a9a9a9"))
colpal = unique((AllIsolates %>% arrange(ST))$colorscheme)


# Bin year collected into 4 year periods
AllIsolates = AllIsolates %>% mutate(YearRange = case_when(Year %in% c("1999", "2000","2001","2002") ~ "1999-2002",
                                                           Year %in% c("2003","2004","2005","2006") ~ "2003-2006",
                                                           Year %in% c("2007", "2008","2009","2010") ~ "2007-2010",
                                                           Year %in% c("2011","2012","2013","2014") ~ "2011-2014",
                                                           Year %in% c("2015", "2016","2017-18") ~"2015-2018", 
                                                           Year=="2022" ~"2022"))

######################################
# Figure 1B -- proportion of each ST
######################################

# Find middle of each year bin 
AllIsolates = AllIsolates %>% mutate(ModifiedYear = case_when(YearRange =="1999-2002" ~ 2000.5, 
                                                              YearRange == "2003-2006" ~ 2004.5,
                                                              YearRange =="2007-2010" ~ 2008.5,
                                                              YearRange=="2011-2014" ~ 2012.5,
                                                              YearRange=="2015-2018" ~ 2016.5,
                                                              YearRange=="2022" ~ 2022,
                                                              YearRange=="2023" ~ 2023))


AllIsolates$countnum=1

# Deduplicate to first observation per patient
FirstObservation = AllIsolates %>% filter(IncludeDeduplicated == 1)

AllIsolatesGrouped = FirstObservation  %>% group_by(YearRange) %>% mutate(numYearRange=sum(countnum)) %>% ungroup() %>% group_by(YearRange, ST) %>% summarize(PropSum=100*(sum(countnum)/mean(numYearRange)), ModifiedYear=ModifiedYear,ST=ST) %>% unique()

AllCombos = FirstObservation %>% select(ST, YearRange) %>% tidyr::expand(ST, YearRange)

AllIsolatesGroupedMod = AllCombos %>% left_join(AllIsolatesGrouped,by=c("ST","YearRange"))
AllIsolatesGroupedMod$PropSum = if_else(is.na(AllIsolatesGroupedMod$PropSum),0, AllIsolatesGroupedMod$PropSum)

AllIsolatesGroupedMod$PropSum = if_else(is.na(AllIsolatesGroupedMod$PropSum),0, AllIsolatesGroupedMod$PropSum)



AllIsolatesGroupedMod = AllIsolatesGroupedMod %>% mutate(ModifiedYear = case_when(YearRange =="1999-2002" ~ 2000.5, 
                                                              YearRange == "2003-2006" ~ 2004.5,
                                                              YearRange =="2007-2010" ~ 2008.5,
                                                              YearRange=="2011-2014" ~ 2012.5,
                                                              YearRange=="2015-2018" ~ 2016.5,
                                                              YearRange=="2022" ~ 2022,
                                                              YearRange=="2023" ~ 2023))
proportionplot <- ggplot(AllIsolatesGroupedMod, aes(x=ModifiedYear, y=PropSum, fill=ST),)+ geom_area(position="fill") +scale_fill_manual(values=colpal)  + scale_x_continuous(breaks=2000:2024) + theme_classic()

# check that everything adds to 100
sum((AllIsolatesGrouped %>% filter(ModifiedYear==2016.5))$PropSum)

######################################
# Figure 1A -- % ST8 and SAE by year
######################################


# % ST8 by year (up to 2022)
#############################
PctDF  =  AllIsolatesGroupedMod %>% filter(ST=="8")
JustPctST8 = ggplot(PctDF, aes(x=ModifiedYear,group=1, y=PropSum)) + geom_point(color=orange,size=3) + geom_line(color=orange,size=1.5) + theme_classic()

#  % SAE per year
#############################

FirstObservationST8s = FirstObservation %>% filter(ST=="8")

# Read in clade assignments based on 
TreeBasedAssignments = read.csv("DataArchive/TreeBasedCladesST8March_V2.csv",header=F)
colnames(TreeBasedAssignments) = c("Identifier","Clade")

IdentifierMap <- read.csv("DataArchive/FullSCL_To_StudyIdentifierMap.csv")
FirstObservationST8s = FirstObservationST8s %>% left_join(IdentifierMap,by="SCL")

FirstObservationST8s = FirstObservationST8s %>% left_join(TreeBasedAssignments,by="Identifier")

FirstObservationSAE = FirstObservationST8s %>% filter(Clade=="SAE")

SAE_by_year = data.frame(YearRange = c(names(table(FirstObservationSAE$YearRange)), "1999-2002", "2003-2006"),NumberSAE = c(table(FirstObservationSAE$YearRange), 0,0))


SAE_by_yearNumTotalYear = FirstObservation %>% group_by(YearRange) %>% summarize(NumYearRange=sum(countnum)) 

SAE_by_year = SAE_by_year %>% left_join(SAE_by_yearNumTotalYear,by="YearRange")


SAE_by_year$PctSAE = (SAE_by_year$NumberSAE/SAE_by_year$NumYearRange)*100
SAE_by_year = SAE_by_year %>% select(YearRange,PctSAE )
SAE_by_year[7,] = c("2023",22.8)

SAE_by_year = SAE_by_year %>% mutate(ModifiedYear = case_when(YearRange =="1999-2002" ~ 2000.5, 
                                                              YearRange == "2003-2006" ~ 2004.5,
                                                              YearRange =="2007-2010" ~ 2008.5,
                                                              YearRange=="2011-2014" ~ 2012.5,
                                                              YearRange=="2015-2018" ~ 2016.5,
                                                              YearRange=="2022" ~ 2022,
                                                              YearRange=="2023" ~ 2023))
SAE_by_year$PctSAE = sapply(SAE_by_year$PctSAE, function(x) as.numeric(as.character(x)))

PctDF_ST8_SAE = SAE_by_year %>% left_join(PctDF %>% select(-ModifiedYear), by="YearRange")


SAE_ST8_yearFixed =  ggplot(PctDF_ST8_SAE, aes(x=ModifiedYear,group=1, y=PropSum)) + geom_point(color=orange,size=3) + geom_line(color=orange,size=1.5) + theme_classic() + geom_line(linetype="dashed", color="purple", size=1, data=PctDF_ST8_SAE, aes(x=ModifiedYear,y=PctSAE)) + geom_point(color="purple", size=3,data=PctDF_ST8_SAE, aes(x=ModifiedYear,y=PctSAE))
SAE_ST8_yearFixed = SAE_ST8_yearFixed + scale_x_continuous(breaks=2000:2023) + guide_legend()


# Pct antibiotic resistant per year
##########################################

PctResistanceDF = AllIsolates %>% filter(IncludeDeduplicated == 1) %>% filter(Genta_int!="")  %>% group_by(YearRange) %>% summarize(NumObs = n())

# Remove those with missing data (just for compositional analyses)
PctResistanceDF


AllIsolatesFirstObs = AllIsolates %>% filter(IncludeDeduplicated == 1)

PctResistanceDF = PctResistanceDF %>%  left_join(AllIsolatesFirstObs %>%filter(Genta_int!="") %>%  group_by(YearRange) %>% filter(FOX_int=="R") %>% summarize(Fox = n()), by="YearRange")
PctResistanceDF = PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>%  group_by(YearRange) %>% filter(CIP_int=="R") %>% summarize(Cip = n()) , by="YearRange")
PctResistanceDF =  PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(Genta_int=="R") %>% summarize(Genta = n()) , by="YearRange")
PctResistanceDF =  PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(Clinda_int=="R") %>% summarize(Clind = n()), by="YearRange")
PctResistanceDF =  PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(Ery_int=="R") %>% summarize(Ery = n()), by="YearRange")
PctResistanceDF =  PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(Rifa_int=="R") %>% summarize(Rifa = n()) , by="YearRange")
PctResistanceDF =  PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(SXT_int=="R") %>% summarize(SXT = n()), by="YearRange")
PctResistanceDF = PctResistanceDF %>% left_join(AllIsolatesFirstObs %>% filter(Genta_int!="") %>% group_by(YearRange) %>% filter(Tetra_int=="R") %>% summarize(Tetra = n()) , by="YearRange")




PctResistanceDF[is.na(PctResistanceDF)] <- 0

PctResistanceDF$Fox = 100*(PctResistanceDF$Fox/PctResistanceDF$NumObs)
PctResistanceDF$Cip = 100*(PctResistanceDF$Cip/PctResistanceDF$NumObs)
PctResistanceDF$Genta = 100*(PctResistanceDF$Genta/PctResistanceDF$NumObs)
PctResistanceDF$Clind = 100*(PctResistanceDF$Clind/PctResistanceDF$NumObs)
PctResistanceDF$Ery = 100*(PctResistanceDF$Ery/PctResistanceDF$NumObs)
PctResistanceDF$Rifa = 100*(PctResistanceDF$Rifa/PctResistanceDF$NumObs)
PctResistanceDF$SXT = 100*(PctResistanceDF$SXT/PctResistanceDF$NumObs)
PctResistanceDF$Tetra = 100*(PctResistanceDF$Tetra/PctResistanceDF$NumObs)

MeltedSusceptibilityYearRange = PctResistanceDF %>% melt(id.vars=c("YearRange", "NumObs" ))
MeltedSusceptibilityYearRange = MeltedSusceptibilityYearRange %>% 
  mutate(ModifiedYear = case_when(YearRange =="1999-2002" ~ 2000.5, 
                                  YearRange == "2003-2006" ~ 2004.5,
                                  YearRange =="2007-2010" ~ 2008.5,
                                  YearRange=="2011-2014" ~ 2012.5,
                                  YearRange=="2015-2018" ~ 2016.5,
                                  YearRange=="2022" ~ 2022,
                                  YearRange=="2023" ~ 2023))

SusceptibilityByYear = ggplot(MeltedSusceptibilityYearRange, aes(x=ModifiedYear, y=value, group=variable, color=variable)) + geom_line(size=1.5) + theme_classic() + scale_color_brewer(palette="Dark2")




# Save plots with shared X axis (year range)
############################################
pdf("FigureOutput/Fig1A-C.pdf",height=12,width=10)
cowplot::plot_grid(SAE_ST8_yearFixed, proportionplot, SusceptibilityByYear+ scale_x_continuous(breaks=2000:2024), align = "v",axis='x',ncol=1)
dev.off()




  
  
  

colorcodes = AllIsolates %>% select(ST, colorscheme) %>% unique()
row.names(colorcodes)= colorcodes$ST


pdf("FigureOutput/S1_ABXSusceptibilityPlots.pdf", width=12,height=8)
for(st in unique(AllIsolatesFirstObs$ST)){
  print(st)
  
  # Go through 'AllIsolates' which has ABX info for all isolates (not just first obs)
  # Make a subset of just ST of interest
  subsetisolates = AllIsolatesFirstObs %>% filter(ST==st)
  subsetisolates[subsetisolates=="" | subsetisolates=="-"| subsetisolates=="_"] <- NA
  subsetisolates = subsetisolates %>% drop_na()
  subsetisolates[subsetisolates=="R"] <- 1
  subsetisolates[subsetisolates=="S" | subsetisolates=="I"] <-0
  
  
  # Make all 1s and 0s numeric  
  subsetisolates[,5:12] <- apply(subsetisolates[,5:12], 2, function(x) as.numeric(as.character(x)))
  
  # Make a 'first row' for cols 5:12 (cefox through tetra) which is the proportion of rows where the col = 1
  firstrow = colSums(subsetisolates[,5:12])/nrow(subsetisolates)
  
  # Row 1 defines the max value
  # Row 2 defines the min value
  # Row 3 is the actual proportion data
  resistanceDF = rbind(rep(1,8), rep(0,8))
  resistanceDF = data.frame(rbind(resistanceDF, firstrow))
  
  colnames(resistanceDF) = c("Cefoxitin","Ciprofloxacin","Gentamicin","Clindamycin","Erythromycin","Rifampicin","Trimethoprim/Sulfamethoxazole","Tetracycline")
  
  STcolor = colorcodes[toString(st),"colorscheme"]
  
  fillcolor = as.vector(c((col2rgb(STcolor)/255)[,1], .75))
  radarchart( resistanceDF[1:ncol(resistanceDF)],pfcol=rgb(fillcolor[1],fillcolor[2],fillcolor[3],fillcolor[4]),caxislabels=c(0,.25,.50,.75,1), title=paste0("ST",st," (n=",toString(nrow(subsetisolates)), ")"),axistype = 1,vlcex=2, calcex=2)
  
  

  }
dev.off()

