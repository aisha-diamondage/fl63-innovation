## tai_plot.py 
This is the streamlit script that performs stAI downstream analysis. 
It allows the user to select significance cut-offs from the RNA-Seq analysis.
counts codon occurrence in the CDS of the genes that are significantly differentially expressed.
Calculates the stairs based on the previously calculated wi's (these were calculated using stAIcalc).
Provides a summary table with the DESeq output, codon count, calculates stAI for the gene's CDS and and the average raw count from all the replicates.
The user is then able to pick any 2 axis from the summary table to use for a scatter plot.
The user can also perform gene enrichment analysis on the differentially expressed genes. 

# The largest heading
## The second largest heading
###### The smallest heading