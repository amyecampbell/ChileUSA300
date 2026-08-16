# Amy Campbell
# Using ggtree to plot SAE subtree 
library(ape) # 5.8
library(ggtree) # v3.10.1
library(ggtreeExtra) # 1.12.0
library(dplyr) # 1.1.4
library(ggstar) # 1.0.4
library(ggplot2) # 3.5.0
library(phytools) # 2.3.0
library(ggnewscale) # 0.4.10

SAEcountries = read.csv("DataArchive/SAEcountries.csv")

##############################################################
# Reference tree (without molecular clock or ancestral recon)
# e.g, not the subtree we'll publish
##############################################################
tree_SAE = ape::read.tree(file="DataArchive/Figure 3/RAxML_bestTree.RaxML_SAE_April_V4_PostCFML")
tree_SAE = ape::root(tree_SAE, outgroup="HUV01",resolve.root=T)

# Make initial circular visualization and extract color palette
TreeVis = ggtree(tree_SAE, layout="circular") + geom_tiplab(offset=.0000005,size=7)
#
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

TreeVis$data = TreeVis$data %>% left_join(SAEcountries, by="label")

colorpal = unique((TreeVis$data %>% arrange(Location,.locale="en"))$color)

# Plot the molecular clock dated tree
#####################################
tree_SAE_bactdating100millionStrict = treeio:::read.beast(file="DataArchive/Figure 3/StrictClock100MillDatesFixed.nex")

MolecClockTree = ggtree(tree_SAE_bactdating100millionStrict, layout="rectangular",right=F) 
MolecClockTreeBackup = ggtree(tree_SAE_bactdating100millionStrict, layout="rectangular",right=F) 

MolecClockTreeRound = ggtree(tree_SAE_bactdating100millionStrict, layout="circular",right=F) 

MolecClockTree$data = MolecClockTree$data %>%  left_join(SAEcountries,by="label")


MolClockTree = ggtree(tree_SAE_bactdating100millionStrict, mrsd="2022-12-02") 
MolClockTree$data = MolClockTree$data %>% left_join(SAEcountries, by="label")
MolClockTree = MolClockTree + theme_tree2() + geom_tippoint(mapping=aes(color=Location)) + scale_color_manual(values=colorpal)



# Collapse nodes containing monophylies entirely of isolates from a single patient 
# Note: MRSA_Chile/ST8/SAETreeVis.R is a private script that contains the code to produce 
# a tree which labels tips by patient ID, enabling the selection of these nodes to
# collapse (monophylies of individual patients' isolates)
nodescollapse = c(262,261,246,242,298,190,198)
MolClockTreeCollapsed=MolClockTree
for(nodeid in nodescollapse){
  MolClockTreeCollapsed = collapse(MolClockTreeCollapsed,node=nodeid)
}

MolClockTreeCollapsed = MolClockTreeCollapsed+ geom_tippoint(mapping=aes(color=Location),size=3) + scale_color_manual(values=colorpal)


# Ancestral reconstruction of geographic states
################################################

# simplest pie plots:
tree_Ancestral = ape::read.tree(file="DataArchive/Figure 3/fullSAEpastMLoutputMolecClock/named.tree_strictClock100MillTree_StrippedDown.nwk")
tree_Ancestral_vis = ggtree(tree_Ancestral, layout="rectangular")
probabilities_Ancestral = read.csv2("DataArchive/Figure 3/fullSAEpastMLoutputMolecClock/marginal_probabilities.character_Location.model_F81.tab",sep="\t")
probabilities_Ancestral[,2:ncol(probabilities_Ancestral)] = apply(probabilities_Ancestral[,2:ncol(probabilities_Ancestral)], 2, function(x) as.numeric(as.character(x)))

pies <- nodepie(probabilities_Ancestral, cols = 2:ncol(probabilities_Ancestral))
pies = pies[tree_Ancestral_vis$data$label]
names(pies) = tree_Ancestral_vis$data$node
pies <- lapply(pies, function(g) g+scale_fill_manual(values = colorpal))


# shows pie probabilities rather than just state with highest marginal posterior prob > 50 
tree_Ancestral_pies <- tree_Ancestral_vis + geom_inset(pies, width=.02, height=.02) 



# Make a mapping of node IDs from the ancestral reconstruction tree to the
# molecular clock tree based on x and y position
###########################################################################
ancestraltreeggobj = ggtree(tree_Ancestral, layout="rectangular",right=F) +geom_tiplab()#+ geom_inset(pies_probs) #+geom_tiplab()

# Make sure the x and y values are to the same # of decimals; x and y positions should be the same between trees, just node labels are diff
ancestraltreeggobj$data$xRounded = round(ancestraltreeggobj$data$x,digits=2)
MolecClockTreeBackup$data$xRounded = round(MolecClockTreeBackup$data$x,digits=2)
ancestraltreeggobj$data$yRounded = round(ancestraltreeggobj$data$y,digits=2)
MolecClockTreeBackup$data$yRounded = round(MolecClockTreeBackup$data$y,digits=2)

# make sure that the x and y values actually match up (they do)
sort(ancestraltreeggobj$data$yRounded) -  sort(MolecClockTreeBackup$data$yRounded)

ancestraltreeggobj$data$xy= paste0(ancestraltreeggobj$data$xRounded, "_", ancestraltreeggobj$data$yRounded)
MolecClockTreeBackup$data$xy= paste0(MolecClockTreeBackup$data$xRounded, "_", MolecClockTreeBackup$data$yRounded)

# make a mapping of node names in the ancestral geographic recon tree and node names in the molecular clock tree (pastML renames them)
mergedTree = data.frame(ancestraltreeggobj$data) %>% left_join(data.frame(MolecClockTreeBackup$data),by="xy")
NodeMappingAncestralToClock = mergedTree %>% select(node.x, node.y, label.x)
colnames(NodeMappingAncestralToClock) = c("node_ancestral","node_clock", "ancestral_label")


# read in a fresh marginal probabilities file 
probabilities_Ancestral = read.csv2("DataArchive/Figure 3/fullSAEpastMLoutputMolecClock/marginal_probabilities.character_Location.model_F81.tab",sep="\t")
probabilities_Ancestral[,2:ncol(probabilities_Ancestral)] = apply(probabilities_Ancestral[,2:ncol(probabilities_Ancestral)], 2, function(x) as.numeric(as.character(x)))
probabilities_Ancestral$Max = apply(probabilities_Ancestral[,2:(ncol(probabilities_Ancestral)-1)], 1,max)
probabilities_Ancestral$ColWithMax = colnames(probabilities_Ancestral[2:ncol(probabilities_Ancestral)])[apply(probabilities_Ancestral[,2:ncol(probabilities_Ancestral)],1,which.max)]

colnames(probabilities_Ancestral) = c("ancestral_label",colnames(probabilities_Ancestral)[2:length(colnames(probabilities_Ancestral))])
probabilities_Ancestral = probabilities_Ancestral %>% left_join(NodeMappingAncestralToClock,by="ancestral_label")

# join so that the molecular clock tree has the mapping of its own nodes to ancestral recon nodes 
MolClockTree$data$node_clock = MolClockTree$data$node
MolClockTree$data = MolClockTree$data %>% left_join(NodeMappingAncestralToClock,by="node_clock")


# Collapse the nodes that are monophyletic clades all from the same patient (label how many were collapsed 
# in a node manually in illustrator)
MolClockTreeWithPiesCollapsed = MolClockTree
MolClockTreePhylo = as.phylo(MolClockTree)
childrennodes = c()
for(nodeid in nodescollapse){
  childrennodes = c(childrennodes, getDescendants(MolClockTreePhylo,nodeid ))
  MolClockTreeWithPiesCollapsed = collapse(MolClockTreeWithPiesCollapsed,node=nodeid)
}

## internal nodes with > 50% marginal probability 
#probabilities_Ancestral_AtLeast50_tips = probabilities_Ancestral %>% filter(Max > .5) %>% filter(node_clock %in% (MolClockTree$data %>% filter(!isTip))$node_clock )

probabilities_Ancestral = probabilities_Ancestral %>% mutate(MaxCountry = if_else(Max>.5, ColWithMax, "NA"))

probabilities_AncestralForJoin = probabilities_Ancestral %>% select(node_clock,MaxCountry )


MolClockTreeWithPiesCollapsed$data = MolClockTreeWithPiesCollapsed$data %>% left_join(probabilities_AncestralForJoin, by="node_clock")


represented= sort(unique((MolClockTreeWithPiesCollapsed$data  %>% filter(!isTip))$MaxCountry))
nodecolorpal = colorkey[represented,"color"]

MolClockTreeWithPiesCollapsed$data
InternalsLabeled = ggtree(MolClockTreeWithPiesCollapsed$data)+ theme_tree2() + geom_nodepoint(aes(color=MaxCountry), shape=15,size=3) +scale_color_manual(values=nodecolorpal)+new_scale_color()+geom_tippoint(aes(color=Location)) + scale_color_manual(values=colorpal)



Figure3_NodesLabeledGeographicState = ggtree(MolClockTreeWithPiesCollapsed$data)+ theme_tree2() + geom_nodepoint(aes(color=MaxCountry), shape=15,size=3) +scale_color_manual(values=nodecolorpal)+new_scale_color()+geom_tippoint(aes(color=Location)) + scale_color_manual(values=colorpal)

# MolClockAncestralRecon_Gr50_nodeslabeled.pdf in original code
ggsave(Figure3_NodesLabeledGeographicState, filename="FigureOutput/Fig3_GeographicMolecularClockTree.pdf",height=8, width=11)

# Chilean circulating clones
##################################
(MolClockTreeWithPiesCollapsed$data %>% filter(node %in% c(214,205,278,235,243,252)))$height_0.95_HPD
(MolClockTreeWithPiesCollapsed$data %>% filter(node %in% c(214,205,278,235,243,252)))$height_0.95_HPD
(MolClockTreeWithPiesCollapsed$data %>% filter(node %in% c(214,205,278,235,243,252)))$node
(MolClockTreeWithPiesCollapsed$data %>% filter(node ==172))$height_0.95_HPD
(MolClockTreeWithPiesCollapsed$data %>% filter(node ==170))$height_0.95_HPD

