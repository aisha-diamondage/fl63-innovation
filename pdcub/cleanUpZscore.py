#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 21:37:03 2021

@author: aisha
"""
import sys
import csv

l=[]
with open("results_up/up_zscores.tsv") as f:
    #reader = csv.reader(f, delimiter='\t')
    for line in f:
        if line not in l:
            print(line.strip())
            l.append(line)