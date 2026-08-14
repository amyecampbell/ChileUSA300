# Amy Campbell
# Allele typing (presence/absence of exact alleles)
# which are diagnostic of
###################

from Bio import SeqIO
import os
import csv

NAEAllelesFiles="/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/NAEAlleles.faa"
SAEAllelesFiles="/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/SAEAlleles.faa"
BothNAESAE="/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/SAE_NAEAlleles.faa"
USA300="/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/USA300Alleles.faa"
PEB="/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/Peb1Alleles.faa"
ProteomesPath="/Users/campbela12/Documents/Planet/ST105/ST8/ST8_Proteomes/FAAs/"


AlleleDict = dict()
KeyList = []
for faa in [NAEAllelesFiles, SAEAllelesFiles, BothNAESAE, USA300,PEB]:

    list_records = list(SeqIO.parse(open(faa),'fasta'))
    for entry in list_records:
        KeyList.append(entry.id)
        AlleleDict[entry.id] = entry.seq

list_proteomes = (os.listdir(ProteomesPath))

OutputArray = [["PPid"]+KeyList]

for proteome in list_proteomes:
    proteomeID=proteome.strip(".faa")
    print(proteomeID)
    sequences = list(SeqIO.parse(open(os.path.join(ProteomesPath, proteome)),'fasta'))
    justseqs=[]
    for s in sequences:
        justseqs.append(str(s.seq))


    smallarray=[proteomeID]

    for k in KeyList:
        if str(AlleleDict[k]) in justseqs:
            smallarray.append("1")
        else:
            smallarray.append("0")
    OutputArray.append(smallarray)
outputfile=open("/Users/campbela12/Documents/Planet/ST105/ST8/BiancoEtAlDiagnostics/AllelePresenceST8.csv","w")
csvwrite = csv.writer(outputfile)

for listitem in OutputArray:
    csvwrite.writerow(listitem)
outputfile.close()
