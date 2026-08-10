#!/bin/bash
# Amy Campbell
# Running QUAST on the Chilean MRSA genomes
source /home/campbela12/miniforge3/bin/activate ST105Env

quast.py -r GCA_008935195.1.fasta -o quastoutput assemblies

