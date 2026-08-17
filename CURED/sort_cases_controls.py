#!/usr/local/bin/python3

import sys 

case_control_file = sys.argv[1]

master_dict = {}
case_list = []
for line in open(case_control_file, 'r'):
    line = line.rstrip().split(',')
    status = line[1]
    genome = line[0].split('.fasta')[0]
    master_dict[genome] = status
    if status == 'case':
        case_list.append(genome)

case_set = set(case_list)
    
report = sys.argv[2]
#print(master_dict)
missing_genomes = []
for line in open(report, 'r'):
    line = line.rstrip().split('\t')
    kmer = line[0]
    controls = []
    local_cases = []
    for item in line[1:]:
        item = item.split(':')[0]
        if master_dict[item] == 'control':
            #acc_list = item.split('_')
            #acc = '_'.join(acc_list[:2])
            item = item + '.fasta'
            controls.append(item)   
        else:
            local_cases.append(item)
    #print(local_cases)
    local_case_set = set(local_cases)
    case_difference = case_set - local_case_set
    #print(len(case_difference))
    #print(case_difference)
    #print(f'{kmer}\t{len(cases)}')
    for G in case_difference:
        if not G in missing_genomes:
            missing_genomes.append(G)
    control_list = '\t'.join(controls)
    print(f'{kmer}\t{control_list}')
    #for c in controls:
     #   c = c + '.fna'
     #   print(c)
#for GENOME in missing_genomes:
#    print(GENOME)
