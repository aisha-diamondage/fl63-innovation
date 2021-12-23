#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 09:41:41 2021

@author: aisha
"""

import sys
import csv

up={}
down={}

with open("neuroAllCAG.tsv") as f:
    reader = csv.reader(f, delimiter = '\t')
    for line in reader:
        if not line[0].startswith("#"):
            regulation = line[0].split("_")[1].strip()
            if regulation == "up":
                if int(line[2].split("-")[0]) <= int(sys.argv[1]):
                    ident = line[1].strip()+"_"+line[2].strip()
                    up[ident] = line
            if regulation == "down":
                if int(line[2].split("-")[0]) <= int(sys.argv[1]):
                    ident = line[1].strip()+"_"+line[2].strip()
                    down[ident] = line
            
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

up_before_count={}
down_before_count={}
up_after_count={}
down_after_count={}

for item in table:
    up_before_count[item] = 0 
    down_before_count[item] = 0
    up_after_count[item] = 0 
    down_after_count[item] = 0
           
            
for chunk in up:
    for neighbor in up[chunk][6].split("|"):
        if neighbor  != ' ':
            aa = neighbor.split(":")[0].split("-")[0].strip()
            if aa == "TAG" or aa == "TAA":
                print(aa, "EQUALS STOP:", up[chunk])
            count = int(neighbor.split(":")[1].split("(")[0].strip())
            up_before_count[aa]= up_before_count[aa] + count
            
for chunk in up:   
    for neighbor in up[chunk][7].split("|"):
        if neighbor  != ' ':
             aa = neighbor.split(":")[0].split("-")[0].strip()
             #if aa == "TAG" or aa == "TAA":
              #  print(aa, "EQUALS STOP:",up[chunk])
             count = int(neighbor.split(":")[1].split("(")[0].strip())
             up_after_count[aa]= up_after_count[aa] + count

#print(up_before_count)
#print(up_after_count)


for chunk in down:
    for neighbor in down[chunk][6].split("|"):
        if neighbor  != ' ':
            aa = neighbor.split(":")[0].split("-")[0].strip()
            if aa == "TAG" or aa == "TAA":
                print(aa, "EQUALS STOP:", down[chunk])
            count = int(neighbor.split(":")[1].split("(")[0].strip())
            down_before_count[aa]= down_before_count[aa] + count

for chunk in down:   
    for neighbor in down[chunk][7].split("|"):
        if neighbor  != ' ':
             aa = neighbor.split(":")[0].split("-")[0].strip()
             #if aa == "TAG" or aa == "TAA":
              #  print(aa, "EQUALS STOP:", down[chunk])
             count = int(neighbor.split(":")[1].split("(")[0].strip())
             down_after_count[aa]= down_after_count[aa] + count


#print(down_before_count)
#print(down_after_count)

CAG_total_up_before=0
CAG_total_up_fater=0
CAG_total_down_before=0
CAG_total_down_after=0

#for item in table:
    
    
 #   print(item, up_before_count[item], up_after_count[item], down_before_count[item], down_after_count[item])