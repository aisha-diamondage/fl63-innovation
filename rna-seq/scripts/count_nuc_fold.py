#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:32:55 2022

@author: aisha
"""

import csv
import glob
import streamlit as st
import os
from Bio import SeqIO
from Bio import Entrez
from Bio import Seq
Entrez.email = 'aisha@diamondage.com'

st.title("RNA-Seq Upset Plots")
st.subheader("This apps shows data from differential analysis of RNA-Seq sample vs vector")

where = st.selectbox("Set", ["set1", "set2"])
pval_cutoff = st.number_input("p-value cutoff", value = 0.05)
logfc_cutoff = st.number_input("log2FC cutoff", value = 1.5)
padj_cutoff = st.number_input("adjusted p-value cutoff", value = 1)
baseMean_cutoff = st.number_input("baseMean cutoff", value = 0)
sort_by = st.radio("Sort by", ["degree", "cardinality"])

print(sort_by)

#where = "set2"
#pval_cutoff = 0.05
#pval_cutoff = 1
#logfc_cutoff = 1.5
#padj_cutoff = 1
#baseMean_cutoff = 0
#sort_by = "degree"
#sort_by = "cardinality"

files = glob.glob(where+"/output/contrasts/*.csv")

d = {}
nucs = []
accessions = []


def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

for file in files:
    #print("NEW LOOP")
    name = os.path.basename(file).split(" ")[0]
    #print(name)
    d[name]=[]
    if name.startswith("Mock") or name.startswith("Vector") or name.startswith("iMet"):
        nucs.append(name)
    else:
        nuc = name[-3:]
        nucs.append(nuc)
    
    with open(file) as f:
        reader = csv.reader(f)
        for i,line in enumerate(reader):
            if i == 0:
                log2fc_index = line.index("log2FC")
                pval_index = line.index("pvalue")
                padj_index = line.index("padj")
                baseMean_index = line.index("baseMean")
                #print(log2fc_index, pval_index,padj_index,baseMean_index )
            else:
                if "NA" in line:
                    continue
                else:
                    if float(line[log2fc_index]) >= logfc_cutoff:
                        if float(line[pval_index]) <= pval_cutoff:
                            if float(line[padj_index]) <= padj_cutoff:
                                if float(line[baseMean_index]) >= baseMean_cutoff:
                                    #print(line)
                                    d[name].append(line[1])
                                    accessions.append(line[1])
                                    #d[line[1]] = line
       
                             
accessions = set(accessions)



### upset plots
from upsetplot import from_contents
from upsetplot import plot
from matplotlib import pyplot

st.set_option('deprecation.showPyplotGlobalUse', False)
upset = from_contents(d)
fig = plot(upset, show_counts=True, orientation= "vertical", sort_by= sort_by )

st.pyplot()





'''

print("fetching fasta records")
#if os.path.exists("seqs.fa"):
 #   os.remove("seqs.fa")
handle = Entrez.efetch(db="nuccore", id= accessions, rettype="fasta_cds_na", retmode="text")
records = SeqIO.parse(handle, "fasta")
for record in records:
    d["_".join(record.id.split("|")[1].split("_")[:2])] = record.description.split(" ")[1].split("=")[1].strip("]")
    
    
    code = [record.seq[i:i+3] for i in range(0, len(record.seq), 3)]
    #print(code)
    for nuc in nucs: 
        try:
            print(nuc, rev_compl(nuc), code.count(rev_compl(nuc)))
        except KeyError:
            print(nuc, 0)
#SeqIO.write(records, "seqs.fa", "fasta")

print (nucs)
'''