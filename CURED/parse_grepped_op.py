#!/usr/local/bin/python3 

import sys 

grepped_op = sys.argv[1]

kmer_dict = {}

for line in open(grepped_op, 'r'):
    line = line.rstrip()
    line_list = line.split()
    kmer = line_list[0].split(':')[1]
    if kmer not in kmer_dict:
        kmer_dict[kmer] = []
    genome_list = line_list[2:]
    for genome in genome_list:
        genome = genome.split(':')[0]
        kmer_dict[kmer].append(genome)

for k,v in kmer_dict.items():
    kmer = k
    genome_set = set(v)
    genomes = '\t'.join(genome_set)
    print(f'{kmer}\t{genomes}')
