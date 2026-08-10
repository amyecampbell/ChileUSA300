# Amy Campbell
# Root ST8 100-rep phylogeny (March V2) for input into TreeShrink to detect
# outlier branches
library(ape) # 5.7.1

ST8Tree = ape::read.tree("/Users/campbela12/Documents/Planet/ST105/ST8/Trees/RAxML_bestTree.RaxML_ST8_March_v2_PostCFML")
ST8Tree = ape::root(ST8Tree, outgroup="USA500", resolve.root=T)
ape::write.tree(ST8Tree, file="/Users/campbela12/Documents/Planet/ST105/ST8/Trees/RAxML_bestTree.RaxML_ST8_March_v2_PostCFML_rootedUSA500.newick")
