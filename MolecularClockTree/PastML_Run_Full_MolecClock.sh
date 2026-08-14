#!/bin/sh
source /home/campbela12/miniforge3/bin/activate pastmlenv
# PastML on the full SAE tree by country 

pastml --tree ../DataArchive/Figure\ 3/strictClock100MillTree_StrippedDown.newick --data ./DataArchive/SAEcountries.csv --columns Location --html_compressed sae_fulltreeMolecClock.html --data_sep , --html sae_fulltree_fullMolecClock.html --work_dir ../DataArchive/Figure\ 3/fullSAEpastMLoutputMolecClock
