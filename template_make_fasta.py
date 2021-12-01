from Bio import SeqIO
import pandas as pd
import textwrap
import os
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




print("creating MP006 nuc fasta")
if os.path.exists("MP006_high_expression_nuc_readthrough.fasta"):
  os.remove("high_expression_nuc_readthrough.fasta")
with open("MP006_nuc_readthrough.fasta", "a") as MP006_nuc_file:
    MP006_nuc_file.write("#>CDS_or_readthrough(RT)|genbankID|description|length_of_nuc_seq")
    for i,record in enumerate(SeqIO.parse("high_expression_gb_seqs.gb", "genbank")):

        for feature in record.features:
            if feature.type == "CDS":
                if str(feature.location).split("(")[1] != "+)":
                    head_error=[">error", record.id, record.description]
                    print("|".join(head_error))#, feature.location)
                    print("This record is has non-standard CDS coordinates:", feature.location)
                else:
                    CDS_start=int(str(feature.location).split(":")[0].strip("["))
                    CDS_end=int(str(feature.location).split(":")[1].split("]")[0])
                    CDS_len=len(str(record.seq)[int(CDS_start):int(CDS_end)])
                    RT_len=len(str(record.seq)[int(CDS_start):])
                    head_CDS=["\n>CDS", record.id, "_".join(str(record.description).split(" ")), str(CDS_len)+'\n']
                    head_RT=["\n>RT", record.id, "_".join(str(record.description).split(" ")), str(RT_len)+'\n']
                    #print("|".join(head_CDS))
                    #print('\n'.join(textwrap.wrap(str(record.seq)[CDS_start:CDS_end])))
                    #print("|".join(head_RT))
                    #print('\n'.join(textwrap.wrap(str(record.seq)[CDS_start:len(str(record.seq))])))

                    MP006_nuc_file.write("|".join(head_CDS))
                    MP006_nuc_file.write('\n'.join(textwrap.wrap(str(record.seq)[CDS_start:CDS_end])))
                    MP006_nuc_file.write("|".join(head_RT))
                    MP006_nuc_file.write('\n'.join(textwrap.wrap(str(record.seq)[CDS_start:len(str(record.seq))])))
MP006_nuc_file.close()

aa=[]
for codon in table:
    aa.append(table[codon])

aa=list(set(aa))
aa.remove("_")
#print(len(aa), aa)

print("creating MP006 prot fasta")
if os.path.exists("MP006_prot_readthrough.fasta"):
  os.remove("MP006_prot_readthrough.fasta")
with open("MP006_prot_readthrough.fasta", "a") as tRNA_prot_file:
    tRNA_prot_file.write("#>CDS_or_readthrough(RT)|genbankID|description|length_of_protein")
    for record in SeqIO.parse("MP006_nuc_readthrough.fasta", "fasta"):
        seq=str(record.seq)
        protein =""
        if (record.id).startswith("CDS"):
            if len(seq)%3 == 0:
                for i in range(0, len(seq), 3):
                    codon = seq[i:i + 3]
                    protein+= table[codon]
                head=str(record.id).split("|")
                head[0]="\n>"+head[0]
                head[3]=str(len(protein.strip("_")))+"\n"
                tRNA_prot_file.write("|".join(head))
                tRNA_prot_file.write('\n'.join(textwrap.wrap(protein.strip("_"))))

        if (record.id).startswith("RT"):
            remainder=len(seq)%3
            seq=seq[:len(seq)-remainder]
            if len(seq)%3 == 0:
                for i in range(0, len(seq), 3):
                    codon = seq[i:i + 3]
                    protein+= table[codon]
                protein = protein.split("_")
                protein = list(filter(None, protein))
                protein = [protein[0], protein[1]]
                protein = 'R'.join(protein)
                #print(protein.split("_"))
                head=str(record.id).split("|")
                head[0]="RT-R"
                head[0]="\n>"+head[0]
                head[3]=str(len(protein))+"\n"
                tRNA_prot_file.write("|".join(head))
                tRNA_prot_file.write('\n'.join(textwrap.wrap(protein)))
print("Done!")
