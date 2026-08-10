#!/bin/bash
#conda activate TreeShrinkEnv

Rscript RootST8Phylogeny.R
run_treeshrink.py -q "0.10" -t /Users/campbela12/Documents/Planet/ST105/ST8/Trees/RAxML_bestTree.RaxML_ST8_March_v2_PostCFML_rootedUSA500.newick -m per-gene -o TreeLengthDetectionST8MarchV2
