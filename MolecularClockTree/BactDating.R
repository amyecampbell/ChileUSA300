# Amy Campbell
# BactDating analysis 
# R 3.4.2

library("BactDating") # 1.1.2
library("ape") # 5.7.1
library("coda") # 0.19.4.1
library("ggtree") # 3.10.0
library("Biostrings") # 2.70.3
library("treeio") # 1.26.0
library("dplyr") # 1.1.4
set.seed(19104)

tree_SAE = ape::read.tree(file="./DataArchive/Figure 2/RAxML_bestTree.RaxML_SAE_April_V4_PostCFML")
tree_SAE = ape::root(tree_SAE, outgroup="HUV01",resolve.root=T)
ggtree(tree_SAE, layout="circular") + geom_tiplab()
dateinfo = read.csv("./DataArchive/Figure 3/SAE_Dated_GenomesDates_Locations.csv")
dateinfo = dateinfo %>% filter(!(is.na(StartRange)))
dateinfo = dateinfo %>% select(Assembly, StartRange, EndRange) %>% unique()
row.names(dateinfo) = dateinfo$Assembly
tree_SAE$edge.length = tree_SAE$edge.length*2825670 # 2825670 sites

dateinfo_bact = cbind(dateinfo$StartRange,dateinfo$EndRange)
row.names(dateinfo_bact) = row.names(dateinfo)

for_root_to_tip = rowMeans(dateinfo_bact)


root_tip=roottotip(tree_SAE,for_root_to_tip)


#Alignment = Biostrings::readDNAMultipleAlignment("/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/SAE_CFML_V4.filtered.fasta", "fasta")
# 2851350 bases in V1, 2825648 in V2, 2825670 in V3, 2825670 in V4


start_arc100 = timestamp()
BactDateRun_arc_dontchangeroot100 = bactdate(tree_SAE,dateinfo_bact,nbIts=1e8, model="arc",showProgress = T, updateRoot = F)
end_arc100 = timestamp()
#save(BactDateRun_arc_dontchangeroot100, file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/arcrun100million_RootFixed_April_V4.rda")

start_strict100 = timestamp()
BactDateRun_strict_dontchangeroot100 = bactdate(tree_SAE,dateinfo_bact,nbIts=1e8, model="strictgamma",showProgress = T, updateRoot = F)
end_strict100 = timestamp()
#save(BactDateRun_strict_dontchangeroot100, file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/strictgammarun100million_RootFixed_April_V4.rda")

start_mixed100 = timestamp()
BactDateRun_mixed_dontchangeroot100 = bactdate(tree_SAE,dateinfo_bact,nbIts=1e8, model="mixedgamma",showProgress = T, updateRoot = F)
end_mixed100 = timestamp()
#save(BactDateRun_strict_dontchangeroot100, file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/mixedgammarun100million_RootFixed_April_V4.rda")


# resetting seed so I can run this at a later date than the previous ones
set.seed(19143)
startrelaxedgamma = timestamp()
BactDateRun_relaxedgamma_dontchangeroot100 = bactdate(tree_SAE,dateinfo_bact,nbIts=1e8, model="relaxedgamma",showProgress = T, updateRoot = F)
endpoiss= timestamp()
save(BactDateRun_relaxedgamma_dontchangeroot100, file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/relaxedgammrun100million_RootFixed_April_V4.rda")

startpoiss = timestamp()
BactDateRun_poiss_dontchangeroot100 = bactdate(tree_SAE,dateinfo_bact,nbIts=1e8, model="poisson",showProgress = T, updateRoot = F)
endpoiss= timestamp()
save(BactDateRun_poiss_dontchangeroot100, file="/Users/campbela12/Documents/Planet/ST105/ST8/MolecularClock/poissonrun100million_RootFixed_April_V4.rda")


# ARC had following ESS values:
# mu: 501
# sigma: 440.2391
# alpha: 501
coda_format_arc100 = as.mcmc.resBactDating(BactDateRun_arc_dontchangeroot100)
coda::effectiveSize(coda_format_arc100_2)
plot(BactDateRun_arc_dontchangeroot100, type="trace")
# Strict gamma had the following ESS values
# mu: 501
# sigma: 0 
# alpha: 424.3712
coda_format_strict100_2 = as.mcmc.resBactDating(BactDateRun_strict_dontchangeroot100)
coda::effectiveSize(coda_format_strict100_2)

# Mixed gamma had the following ESS values
# mu: 501
# sigma: 501  
# alpha: 501
coda_format_mixed100_2 = as.mcmc.resBactDating(BactDateRun_mixed_dontchangeroot100)
coda::effectiveSize(coda_format_mixed100_2)

# Relaxed gamma had hte following ESS values:
# 501, 501, 501
coda_format_relaxed100_2 = as.mcmc.resBactDating(BactDateRun_relaxedgamma_dontchangeroot100)
coda::effectiveSize(coda_format_relaxed100_2)

# Strict (poisson, so discrete) had the following ESS values:
# mu: 595.0917
# sigma: 0
# alpha: 1371.9059
coda_format_strict100_2 = as.mcmc.resBactDating(BactDateRun_poiss_dontchangeroot100,burnin = .5)
coda::effectiveSize(coda_format_strict100_2)


modelcompare(BactDateRun_mixed_dontchangeroot100,BactDateRun_arc_dontchangeroot100)
modelcompare(BactDateRun_arc_dontchangeroot100,BactDateRun_strict_dontchangeroot100)
modelcompare(BactDateRun_strict_dontchangeroot100,BactDateRun_mixed_dontchangeroot100)
modelcompare(BactDateRun_strict_dontchangeroot100,BactDateRun_relaxedgamma_dontchangeroot100)
modelcompare(BactDateRun_poiss_dontchangeroot100,BactDateRun_relaxedgamma_dontchangeroot100)
modelcompare(BactDateRun_poiss_dontchangeroot100,BactDateRun_strict_dontchangeroot100)


# Extracting strict gamma molec clock tree for better visualization
###########################################
plot(BactDateRun_strict_dontchangeroot100, type="treeCI")


highestdate = max(dateinfo$EndRange)
treedata_strict = as.treedata.resBactDating(BactDateRun_strict_dontchangeroot100)
treedata_strict_extracted = methods::new('treedata',phylo=treedata_strict[[1]],data=dplyr::tbl_df(as.data.frame(treedata_strict[[2]])))
treedata_strict_extracted@data$height_0.95_HPD = lapply(treedata_strict_extracted@data$height_0.95_HPD, function(x) highestdate-x)
write.beast(treedata_strict_extracted,'DataArchive/Figure 3/StrictClock100MillDatesFixed.nex')

# Subsequently edited in FigTree v1.4.4 
write.beast.newick(treedata_strict_extracted,'DataArchive/Figure 3/StrictClock100MillDatesFixed.newick')


# 