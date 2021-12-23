#to run, place in folder with CDS fasta or with csv containing NM accessions and one header line.
#python count_CAG.py [#int_segment_size] [#int_sliding_window_size] [#int_%significant]

from Bio import SeqIO
import pandas as pd
import textwrap
import os
import sys
from Bio import Entrez
Entrez.email = 'aisha@diamondage.com'
import glob
import csv

#### CHANGE PARAMETERS TO TEST ####
test = int(sys.argv[1])
slide = int(sys.argv[2])
percent_sig = int(sys.argv[3])


files = glob.glob("*csv")


#print("checking fasta files")
for file in files:
    file_name = file.strip(".csv")
    #print(file_name)
    genes=[]

    if os.path.exists(file_name+".fasta"):
        #print("record exists.")
        continue
    else:
        #print("fetching record.")
        with open(file) as f:
            reader = csv.reader(f)
            for i, line in enumerate(reader):
                if i > 0:
                    genes.append(line[0])

            handle = Entrez.efetch(db="nuccore", id=genes, rettype="fasta_cds_na", retmode="text")
            #print(handle)
            #Entrez.efetch(db="nuccore", id= high_expression_genes, rettype="gb, retmode="text")
            records = SeqIO.parse(handle, "fasta")
            print(records)
            SeqIO.write(records, file_name+".fasta", "fasta")
            handle.close()
#print("Done!")

table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


print("#file","\t","gene","\t","chunk_coordinates","\t","chunk_length","\t","CAG_count","\t","CAG%","\t","before_neighbours","\t","after_neighbours")
files = glob.glob("Neur*.fasta")
for file in files:
    file_name = file.strip(".fasta")
    #print("genes in", file)
    for record in SeqIO.parse(file, "fasta"):
        gene = record.name.split("|")[1].split("_cds_")[0]
        #print(gene)
        seq_codons=[]
        seq=record.seq
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                seq_codons.append(str(codon))
                #print(seq_codons)
            for i in range(0, len(seq_codons)-test+slide, slide):
                before=[]
                after=[]
                chunk = seq_codons[i:i + test]
                chunk_len = len(chunk)
                bin = str(i)+"-"+str(i+chunk_len)
                cag_count = chunk.count("CAG")
                cag_percent = round(100*(chunk.count("CAG")/len(chunk)),2)
                if gene == " NM_144734.1":
                    print(gene, chunk)
                #extract neighbours, one before and one after_count
                #unless the CAG is on the edge of the chunk, it will have one neighbour
                for i,codon in enumerate(chunk):
                    #print(i, chunk)
                    #print(i, codon, len(chunk))
                    if codon == "CAG":
                        if i == 0:
                            #print(chunk, i)
                            #print(chunk[i+1])
                            
                            after.append(chunk[i+1])
                            #print(i, codon, "NA", chunk[i], chunk[i+1])
                        if i == len(chunk)-1:
                            #print(i, chunk)
                            #print(chunk[i-1])
                            before.append(chunk[i-1])
                            #print(i, codon, chunk[i-1], chunk[i], "NA")
                        if i != len(chunk)-1 and i != 0 :
                            #print(i, chunk[i-1])
                            before.append(chunk[i-1])
                            after.append(chunk[i+1])
                            #print(i, codon, chunk[i-1], chunk[i], chunk[i+1])




                #print(i)
                #print(len(seq_codons), len(seq_codons)-slide, len(seq_codons)-test)
                #print("len",chunk_len)
                #print("____________")



                #print only things with significant % CAG
                if cag_percent >= percent_sig:
                    #set up neighbour output
                    before_count=[]
                    after_count=[]
                    for item in table:
                        if before.count(item) > 0:
                            before_count.append(item + "-"+ table[item] + ":" + str(before.count(item)) + "("+ str(round(before.count(item)/len(before)*100, 2)) +"%)")
                        if after.count(item) > 0:
                            after_count.append(item + "-"+ table[item] + ":" + str(after.count(item)) + "("+ str(round(after.count(item)/len(after)*100, 2)) +"%)")

                    #print output
                    #print(len(seq_codons), len(seq_codons)-slide)
                    print(file_name,"\t", gene ,"\t",  bin , "\t", chunk_len, "\t",  cag_count ,"\t",  cag_percent,"\t", "|".join(before_count),"\t", "|".join(after_count))


        else:
            
            print("sequence not divisible by 3, make sure you have the CDS, not gene sequence.")
