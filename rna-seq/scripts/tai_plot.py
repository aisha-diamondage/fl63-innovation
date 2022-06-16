#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed June 15 16:34:55 2022 

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
import pandas as pd
from statistics import geometric_mean

st.title("stAI downstream analysis")
st.subheader("This apps plots the correlation between stAI and codon counts in differentially expressed genes")
             


what = st.selectbox("Sample", ["AsnGTT-3", "GlyTCC-1", "TyrGTA-5", "LeuTAA-1", "ValAAC-2", "ValTAC-1", "ArgTCT-1", "AspGTC-1"])
pval_cutoff = st.number_input("p-value cutoff", value = 0.05)
logfc_cutoff = st.number_input("log2FC cutoff", value = 1.5)
pos_logfc = st.checkbox("Positive log2FC only", value = True)
padj_cutoff = st.number_input("adjusted p-value cutoff", value = 1.0)
baseMean_cutoff = st.number_input("baseMean cutoff", value = 10)


#what = "LeuTAA-1"
#pval_cutoff = 0.05
#pval_cutoff = 1
#logfc_cutoff = 1.5
#padj_cutoff = 1
#baseMean_cutoff = 10
#pos_logfc = True

#cloud dev
#files = glob.glob("tai_data/DE/*.csv")
#tai_file = "tai_data/DE/*.csv"

#local dev
files = glob.glob("../../tai_data/DE/*.csv")
tai_file = "../../tai_data/stAIcalc_out/tRNA_expression_"+what+"/output_wi_file.txt"


d = {}
nucs = []
accessions = []
deseq =  {"accession":["symbol", "codonCount", "log2FC","pvalue","padj","statistic","baseMean", "stAI"]}



############################
# GET SIGNIFICANT DE GENES #
############################


def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

for file in files:
    #print("NEW LOOP")
    name = os.path.basename(file).split(" ")[0]
    if name == what:
        d[name]=[]
    #if name.startswith("Mock") or name.startswith("Vector") or name.startswith("iMet"):
     #   nucs.append(name)
    #else:
     #   nuc = name[-3:]
      #  nucs.append(nuc)
    
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
                        if pos_logfc:
                            if float(line[log2fc_index]) >= logfc_cutoff:
                                if float(line[pval_index]) <= pval_cutoff:
                                    if float(line[padj_index]) <= padj_cutoff:
                                        if float(line[baseMean_index]) >= baseMean_cutoff:
                                            d[name].append(line[1])
                                            accessions.append(line[1])
                                            deseq[line[1]] = line
                        else:
                            if abs(float(line[log2fc_index])) >= logfc_cutoff:
                                if float(line[pval_index]) <= pval_cutoff:
                                    if float(line[padj_index]) <= padj_cutoff:
                                        if float(line[baseMean_index]) >= baseMean_cutoff:
                                            #print(line)
                                            d[name].append(line[1])
                                            accessions.append(line[1])
                                            deseq[line[1]] = line
       
                             
accessions = set(accessions)


#######################
# READ stAIcalc FILES #
#######################

wi = {}
with open(tai_file) as f:
    reader = csv.reader(f, delimiter = '\t')
    for i, line in enumerate(reader):
        if i != 0:
            #print(line[0])
            wi[line[0]] = line[1] 


#######################################################
# GET CODON COUNTS AND GENE SYMBOLS AND CREATE TABLES #
#######################################################

st.subheader("Table of all genes that pass this filter in "+ what + ":")

#counts = {'sample':['accession'], 'gene_name':['rev_codon'] }
#counts = []


#print("*********")
genes = []
handle = Entrez.efetch(db="nuccore", id= accessions, rettype="fasta_cds_na", retmode="text")
records = SeqIO.parse(handle, "fasta")
for record in records:
    gene_acc = ["_".join(record.id.split("|")[1].split("_")[:2])][0]
    gene_name = record.description.split(" ")[1].split("=")[1].strip("]")
    code = [record.seq[i:i+3] for i in range(0, len(record.seq), 3)]
    stai_wi = []
    #counts['sample'].append(gene_acc)
    #counts['gene_name'].append(gene_name)


 
    for sample in d:
        sample_strip = sample.split("-")[0]
        nuc = sample_strip[-3:]
        deseq[gene_acc][0] = gene_name
        deseq[gene_acc][1] = str(code.count(rev_compl(nuc)))
        for trip in code:
            stai_wi.append(float(wi[trip]))
        deseq[gene_acc].append(geometric_mean(stai_wi))
        genes.append(gene_name)
        #if gene_acc in list_to_test_tai_calc:
         #   print(gene_acc, geometric_mean(stai_wi))
        
                
            
            

#transpose = st.checkbox("Transpose")
#df = pd.DataFrame.from_dict(counts).astype(str)
#df = df.set_index('sample')

#print("")
#print(df)

df = pd.DataFrame.from_dict(deseq).astype(str)
df = df.transpose()
df.columns = df.iloc[0]
df = df.drop(df.index[0])

#print(df)

@st.cache
def convert_df(df):
   return df.to_csv().encode('utf-8')


csv = convert_df(df)

st.download_button(
   "Press to Download",
   csv,
   what+".csv",
   "text/csv",
   key='download-csv'
)

    
st.dataframe(df)


#################
# SCATTER PLOTS #
#################


st.subheader("Plots of all genes that pass the above set filters in "+ what + ":")


x_select = st.selectbox("X-axis", ["codonCount", "log2FC","pvalue","padj","statistic","baseMean", "stAI"])
y_select = st.selectbox("Y-axis", ["stAI", "codonCount", "log2FC","pvalue","padj","statistic","baseMean", ])
gene_label = st.checkbox("Gene label", value = True)


#print(c)
#print(df["log2FC"])

x = [float(item) for item in df[x_select].tolist()]
y = [float(item) for item in df[y_select].tolist()]
import matplotlib.pyplot as plt

st.set_option('deprecation.showPyplotGlobalUse', False)
fig = plt.scatter(x,y)
plt.xlabel(x_select)
plt.ylabel(y_select)
if gene_label:
    for i, label in enumerate(df["symbol"].tolist()):
        plt.annotate(label, (x[i], y[i]))
plt.show()

st.pyplot()




############################
# GENE ENRICHMENT ANALYSIS #
############################



st.subheader("Gene enrichment analysis of all genes that pass the above set filters in "+ what + ":")



import gseapy
from gseapy.plot import barplot
#print( gseapy.get_library_name())

GEA_select = st.selectbox("Select gene enrichment analysis database", ['ClinVar_2019' \
 ,'Disease_Perturbations_from_GEO_down' \
 ,'Disease_Perturbations_from_GEO_up' \
 ,'Disease_Signatures_from_GEO_down_2014' \
 ,'Disease_Signatures_from_GEO_up_2014' \
 ,'Elsevier_Pathway_Collection' \
 ,'Gene_Perturbations_from_GEO_down' \
 ,'Gene_Perturbations_from_GEO_up' \
 ,'GO_Biological_Process_2021' \
 ,'GO_Cellular_Component_2021' \
 ,'GO_Molecular_Function_2021' \
 ,'GTEx_Aging_Signatures_2021' \
 ,'GTEx_Tissue_Expression_Down' \
 ,'GTEx_Tissue_Expression_Up' \
 ,'KEGG_2019_Mouse' \
 ,'KEGG_2021_Human' \
 ,'OMIM_Disease' \
 ,'OMIM_Expanded'])

#gseapy.enrichr(gene_list=genes, description='pathway', gene_sets='KEGG_2016', outdir='test')

enr = gseapy.enrichr(gene_list=genes, gene_sets=GEA_select, no_plot=True)

gseapy.plot.barplot(enr.results, column='Adjusted P-value', title='', ofname="gea")
        
csv = convert_df(enr.res2d)
st.download_button(
   "Press to Download",
   csv,
   what+"_"+GEA_select+"_GEA.csv",
   "text/csv",
   key='download_GEA'
)

enr.res2d

st.image("gea.png")

