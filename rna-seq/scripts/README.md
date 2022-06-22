
## count_nuc_fold.py

### This is a streamlit app that is used to create upset plots.
This tool is menat to answer the question: are the genes that are differentially expressed between "treatment" and vector unique to each treatment, or are they artifacts and appear in most differentiall analysis? 


- It allows the user to select significance cut-offs from the RNA-Seq analysis (treatment v vector).
- It creates upset plots of all genes that pass the set cut-offs in ALL treatment conditions.
- Counts all codon occurrence in the CDS of the genes for all treatment codons in the genes that are differentially expressed. 
- Provides a table summary of ALL genes that are significantly expressed across all treatment conditions vs vector. In this table, the count of each codon in each gene is provided, along with a T/F notations for whether this gene was significantlly differentially expressed int his treatment condition. 




## tai_plot.py 

### This is the streamlit script that performs stAI downstream analysis. 
This tool is meant to answer the question: is there any correlation or enrichment in the genes that are differentially expressed between our "treatment" conditions and the vector?


- It allows the user to select a treatment condition.
- The user then selects significance cut-offs from the RNA-Seq analysis (treatment v vector).
- Counts codon occurrence in the CDS of the genes that are significantly differentially expressed.
- Calculates the stairs based on the previously calculated wi's (these were calculated using stAIcalc).
- Provides a summary table with the DESeq output, codon count, calculates stAI for the gene's CDS and and the average raw count from all the replicates.
- The user is then able to pick any 2 axis from the summary table to use for a scatter plot.
- The user can also perform gene enrichment analysis on the differentially expressed genes. 


## requirements.txt

This file lists the requirements needed to launch the streamlit app. It's for both the tai_plot app and count_nuc_fold app.

