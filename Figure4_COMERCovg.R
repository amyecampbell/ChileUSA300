# Amy Campbell
# Figure 4A: Tree of SAE genomes with 
# finding % covered by HSPS with >99% id in COMER
library(dplyr) # 1.1.4
library(stringr) # 1.5.1
library(ggplot2) # 3.5.0
library(reshape2) # 1.4.4
library(ggtree) # 3.10.1
library(ggnewscale) # 0.4.10
library(phytools) # 2.3.0
library(ggtreeExtra) # 1.12.0

##########################
# Blast results for COMER
##########################

colnamesblast = c("name","qseqid","slen","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs","qcovhsp")

COMER_Blast = read.csv2("./DataArchive/COMERall.tab", sep="\t",header=F)

TreeBasedClades= read.csv("./DataArchive/TreeBasedCladesST8March_V2.csv",header=F)
colnames(TreeBasedClades) = c("Genome","Clade")

TreeBasedClades$Genome = sapply(TreeBasedClades$Genome, function(x) str_replace(x, " ","_"))


# % COMER covered with any combination of HSPs having 99% identity covered by genome
colnames(COMER_Blast) = colnamesblast

COMER_Blast$pident = sapply(COMER_Blast$pident,function(x) as.numeric(as.character(x)))
COMER_Blast$length = sapply(COMER_Blast$length,function(x) as.numeric(as.character(x)))
COMER_Blast$qstart = sapply(COMER_Blast$qstart,function(x) as.numeric(as.character(x)))
COMER_Blast = COMER_Blast %>% filter(pident>99)

AllCols = unique(COMER_Blast$name)
BlastResultsCOMER = data.frame(matrix(0, 24688, length(AllCols)))
colnames(BlastResultsCOMER) = AllCols


for(row in 1:nrow(COMER_Blast)){
  startbp = (COMER_Blast[row, "qstart"])
  endbp = (COMER_Blast[row, "qend"])
  seqid =  (COMER_Blast[row, "name"])
  #print(seqid)
  BlastResultsCOMER[startbp:endbp,c(seqid)] <- 1
}
PctCoveredCOMER = 100*(colSums(BlastResultsCOMER)/24688)

BlastResultsCOMER$position=1:nrow(BlastResultsCOMER)

Missing = setdiff(TreeBasedClades$Genome , colnames(BlastResultsCOMER))
MissingMatrix = data.frame(matrix(0, 24688, length(Missing)))
colnames(MissingMatrix) = Missing

BlastResultsCOMERMerged = cbind(BlastResultsCOMER,MissingMatrix )
ForHeatMap = BlastResultsCOMERMerged %>% melt(id.vars=c("position"))

colnames(ForHeatMap) = c("Position","Genome", "Covered99")
ordered = TreeBasedClades %>% arrange(Clade)


# Pct of COMER covered in just SAE genomes
##########################################
Included_SAE_Tree = read.csv("./DataArchive/SAE_Dated_GenomesDates_Locations.csv")
PP_SCL_Map = read.csv("./DataArchive/FullSCL_To_StudyIdentifierMap.csv")

JustSAEs = ForHeatMap %>% filter(Genome %in% Included_SAE_Tree$Assembly )

tree_SAE = ape::read.tree(file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/RAxML_bestTree.RaxML_SAE_April_V4_PostCFML")
tree_SAE = ape::root(tree_SAE, outgroup="MRSA_S9",resolve.root=T)
tree_SAE= ape::drop.tip(tree_SAE,tip = "Reference")
OrderTips = (ggtree(tree_SAE))$data %>% arrange(y) %>% filter(isTip)

heatmapCOMERsae=ggplot(JustSAEs, aes(x=Position, y=Genome, fill=factor(Covered99))) + geom_tile() +theme(axis.text.y=element_text(size=3,vjust=.5))+ scale_fill_manual(values=c("white","black"))

setdiff(OrderTips$label, Included_SAE_Tree$Assembly)

heatmapCOMERsae$data$Genome = factor(heatmapCOMERsae$data$Genome, levels=unique(OrderTips$label))
ggsave(heatmapCOMERsae, file="FigureOutput/Figure4_CovgCOMERsae.png",width=20,height=15)


##############################################################
# Reference tree (without molecular clock or ancestral recon)
##############################################################
SAEcountries = read.csv("./DataArchive/SAEcountries.csv")

tree_SAE = ape::read.tree(file="./DataArchive/RAxML_bestTree.RaxML_SAE_April_V4_PostCFML")


tree_SAE = ape::root(tree_SAE, outgroup="HUV01",resolve.root=T)




TreeVis = ggtree(tree_SAE, layout="circular") + geom_tiplab(offset=.0000005,size=7)


SAEcountries = SAEcountries %>% mutate(color=case_when(Location=="Chile" ~"#DDCC77",
                                                       Location=="Bolivia" ~ "#AA4499",
                                                       Location=="Colombia" ~"#332288",
                                                       Location=="Argentina" ~ "#CC6677",
                                                       Location=="US" ~ "#44AA99",
                                                       Location=="Germany"~"#6699CC",
                                                       Location=="UK" ~ "#88CCEE",
                                                       Location=="Australia"~ "#117733",
                                                       Location=="Uruguay"~ "#661100",
                                                       Location=="Paraguay" ~"#CAB0E2",
                                                       Location=="Ecuador"~ "#999933"))

colorkey = SAEcountries %>% select(color, Location) %>% unique()
rownames(colorkey) = colorkey$Location
SAEcountries$label = SAEcountries$Assembly
TreeVisRect = ggtree(tree_SAE, layout="rectangular", branch.length="none") 
TreeVisRect2 = ggtree(tree_SAE, layout="rectangular") 

TreeVisRect$data = TreeVisRect$data %>% left_join(SAEcountries,by="label")

TreeVisRect = TreeVisRect + geom_tiplab(size=8) +geom_tippoint(mapping=aes(color=Location), 
                                                               size=2.8,
                                                               show.legend=T) + scale_color_manual(values=colorpal)

ggsave(TreeVisRect, file="./FigureOutput/Figure4_RectangularTree.pdf",height=15, width=10)

