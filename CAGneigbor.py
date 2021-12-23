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
        if neighbor.strip()  != '':
            aa = neighbor.split(":")[0].split("-")[0].strip()
            if aa == "TAG" or aa == "TAA":
                print(aa, "EQUALS STOP:", up[chunk])
            #print(neighbor.split(":"))
            count = int(neighbor.split(":")[1].split("(")[0].strip())
            up_before_count[aa]= up_before_count[aa] + count
            
for chunk in up:   
    for neighbor in up[chunk][7].split("|"):
        if neighbor.strip()   != '':
             aa = neighbor.split(":")[0].split("-")[0].strip()
             #if aa == "TAG" or aa == "TAA":
              #  print(aa, "EQUALS STOP:",up[chunk])
             count = int(neighbor.split(":")[1].split("(")[0].strip())
             up_after_count[aa]= up_after_count[aa] + count

#print(up_before_count)
#print(up_after_count)


for chunk in down:
    for neighbor in down[chunk][6].split("|"):
        if neighbor.strip()   != '':
            aa = neighbor.split(":")[0].split("-")[0].strip()
            if aa == "TAG" or aa == "TAA":
                print(aa, "EQUALS STOP:", down[chunk])
            count = int(neighbor.split(":")[1].split("(")[0].strip())
            down_before_count[aa]= down_before_count[aa] + count

for chunk in down:   
    for neighbor in down[chunk][7].split("|"):
        if neighbor.strip()   != '':
             aa = neighbor.split(":")[0].split("-")[0].strip()
             #if aa == "TAG" or aa == "TAA":
              #  print(aa, "EQUALS STOP:", down[chunk])
             count = int(neighbor.split(":")[1].split("(")[0].strip())
             down_after_count[aa]= down_after_count[aa] + count


#print(down_before_count)
#print(down_after_count)

CAG_total_up_before=0
CAG_total_up_after=0
CAG_total_down_before=0
CAG_total_down_after=0

for item in table:
    CAG_total_up_before = CAG_total_up_before + up_before_count[item]
    CAG_total_up_after = CAG_total_up_after + up_after_count[item]
    CAG_total_down_before = CAG_total_down_before + down_before_count[item]
    CAG_total_down_after = CAG_total_down_after + down_after_count[item]

#print(CAG_total_up_before, CAG_total_up_after, CAG_total_down_before,CAG_total_down_after )

# tAI scores for codons in humans (scores retrieved from Tuller 2010 translational ramp paper)
tAI = {}
tAI["TTT"] = 0.161002
tAI["TTC"] = 0.366748
tAI["TTA"] = 0.213936
tAI["TTG"] = 0.282396
tAI["TCT"] = 0.336186
tAI["TCC"] = 0.242054
tAI["TCA"] = 0.152845
tAI["TCG"] = 0.171149
tAI["TAT"] = 0.218399
tAI["TAC"] = 0.449878
tAI["TAA"] = 0.061128
tAI["TAG"] = 0.050122
tAI["TGT"] = 0.402506
tAI["TGC"] = 0.91687
tAI["TGA"] = 0.091687
tAI["TGG"] = 0.304401
tAI["CTT"] = 0.366748
tAI["CTC"] = 0.264059
tAI["CTA"] = 0.091724
tAI["CTG"] = 0.334963
tAI["CCT"] = 0.305623
tAI["CCC"] = 0.220049
tAI["CCA"] = 0.213967
tAI["CCG"] = 0.190709
tAI["CAT"] = 0.147586
tAI["CAC"] = 0.336186
tAI["CAA"] = 0.336186
tAI["CAG"] = 0.749389
tAI["CGT"] = 0.213936
tAI["CGC"] = 0.154034
tAI["CGA"] = 0.183395
tAI["CGG"] = 0.211491
tAI["ATT"] = 0.535208
tAI["ATC"] = 0.552567
tAI["ATA"] = 0.152855
tAI["ATG"] = 0.611247
tAI["ACT"] = 0.305623
tAI["ACC"] = 0.220049
tAI["ACA"] = 0.183405
tAI["ACG"] = 0.242054
tAI["AAT"] = 0.459902
tAI["AAC"] = 1
tAI["AAA"] = 0.519563
tAI["AAG"] = 0.685819
tAI["AGT"] = 0.107335
tAI["AGC"] = 0.244499
tAI["AGA"] = 0.183374
tAI["AGG"] = 0.211491
tAI["GTT"] = 0.336186
tAI["GTC"] = 0.242054
tAI["GTA"] = 0.152845
tAI["GTG"] = 0.537897
tAI["GCT"] = 0.886308
tAI["GCC"] = 0.638142
tAI["GCA"] = 0.27515
tAI["GCG"] = 0.240831
tAI["GAT"] = 0.254921
tAI["GAC"] = 0.580685
tAI["GAA"] = 0.397311
tAI["GAG"] = 0.52445
tAI["GGT"] = 0.201253
tAI["GGC"] = 0.458435
tAI["GGA"] = 0.275061
tAI["GGG"] = 0.301956

print("#The tRNA adaptation index (tAI) is a widely used measure of the efficiency by which a coding sequence is recognized by the intra-cellular tRNA pool")
print("# tAI scores for codons in humans were retreived from the PDCUB BioArchives paper in which they use scores retrieved from Tuller 2010 translational ramp paper\n")
print("#codon", "\t", "tAI_score", "\t", "up_before","\t", "up_after","\t",  "down_before","\t",  "down_after")

for item in table:   
    ub = round(100*(up_before_count[item]/CAG_total_up_before), 2)
    ua = round(100*(up_after_count[item]/CAG_total_up_after), 2)
    db = round(100*(down_before_count[item]/CAG_total_down_before), 2)
    da = round(100*(down_after_count[item]/CAG_total_down_after), 2)
    tai_score= round(tAI[item], 2)
    
    print(item,"\t",tai_score , "\t", ub,"\t", ua ,"\t", db ,"\t", da )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    