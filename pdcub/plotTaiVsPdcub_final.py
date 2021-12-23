import sys,re,os,getopt
import numpy as np
from numpy import mean
from numpy import array
from Bio import SeqIO
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('Agg') # enables saving graphics
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict
import math

###############
# SUBROUTINES #
###############

# PDCUB scores file looks like:
# ENST00000641515.2       OR4F5   60      59      0       2618    2.9214349680077714
def readPdcubScores(pdcubScores):
    pdcubScoresDict = {}
    f1 = open(pdcubScores)
    dataf1 = f1.readlines()
    linePattern = re.compile('(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)')
    for line in dataf1:
        if linePattern.search(line):
            match = linePattern.search(line)
            shortID = match.group(1)
            pdcubScore = match.group(7)
            pdcubScoresDict[shortID] = pdcubScore
            #print(dcubScoresDict)
    f1.close()
    return pdcubScoresDict

def validTranscript(cdsLength,cds,initiator,stop):
    if cdsLength % 3 == 0 and cdsLength >= (binCountNt*binSizeNt) and initiator == "ATG" and stop in stops:
        return True
    else:
        return False

def processFastaFile(fastaFile):
    #cdsPattern = re.compile('CDS:(\d+)\-(\d+)')
    #geneRecordPattern = re.compile('(ENST\d*\.\d*)\S*')
    sequences = SeqIO.parse(fastaFile,"fasta")
    valid_mRNAs = []
    for record in sequences:
        
        description = record.description
        transcript = record.seq
        transcript = transcript.upper()
        #transcript = transcript.transcribe()
        geneID = str(record.id).split('|')[1]
        # geneID looks like: ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094|OTTHUMT00000003223.1|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
        #shortMatch = geneRecordPattern.search(geneID)
        shortID = str(record.id).split('|')[1].split("_cds_")[0]
        #if cdsPattern.search(description):
           # match = cdsPattern.search(description)
        start = 0
        end = len(record.seq)
        cds = record.seq
        cdsLength = len(cds)
        initiator = transcript[start:start+3]
        stop = transcript[end-3:end+1]
        #print(initiator, stop)
        if validTranscript(cdsLength,cds,initiator,stop):
                #print(cdsLength,cds,initiator,stop)
                if shortID in pdcubScoresDict.keys():
                    valid_mRNAs.append((shortID,geneID,str(transcript),str(cds),start,end,cdsLength))
    #print("*VALID MRNA*", valid_mRNAs)
    return valid_mRNAs


def calculateTAIs(valid_mRNAs):

    firstSegmentTAIs = {}
    secondSegmentTAIs = {}
    totalTAIs = {}
    transcriptCumulativePositionalCodons = {}
    cumulativePositionalTAIs = {}
    cumulativeAveragePositionalTAIs = {}
    maxCdsLength = 0

    for mRNA in valid_mRNAs:
        firstSegmentCodons = []
        secondSegmentCodons = []
        totalCodons = []
        shortID,geneID,transcript,cds,start,end,cdsLength = mRNA
        if cdsLength > maxCdsLength:
            maxCdsLength = cdsLength
        transcriptCumulativePositionalCodons[shortID] = {}
        for i in range(0,binSizeNt,3):
            codon = cds[i:i+3]
            firstSegmentCodons.append(str(codon))
        for i in range(binSizeNt,2*binSizeNt,3):
            codon = cds[i:i+3]
            secondSegmentCodons.append(str(codon))
        for i in range(0,2*binSizeNt,3):
            codon = cds[i:i+3]
            totalCodons.append(str(codon))

        # first partial segment of transcript
        codonTai1values = []
        for codonName in firstSegmentCodons:
            codonTai1values.append(tAI[codonName])
        firstSegmentCodonsProduct = float(np.prod(codonTai1values))
        firstSegmentCodonLength = float(len(firstSegmentCodons))
        firstSegmentTAI = firstSegmentCodonsProduct**(1/firstSegmentCodonLength)
        firstSegmentTAIs[shortID] = firstSegmentTAI

        # firstSegmentTAIs looks like this:
        # {'ENST00000494801.5': 0.40805956158458556, 'ENST00000216962.9': 0.3645331933949188, ....

        # second partial segment of transcript
        codonTai2values = []
        for codonName in secondSegmentCodons:
            codonTai2values.append(tAI[codonName])
        secondSegmentCodonsProduct = float(np.prod(codonTai2values))
        secondSegmentCodonLength = float(len(secondSegmentCodons))
        secondSegmentTAI = secondSegmentCodonsProduct**(1/secondSegmentCodonLength)
        secondSegmentTAIs[shortID] = secondSegmentTAI

        # full length
        codonTaiTotalValues = []
        for codonName in totalCodons:
            codonTaiTotalValues.append(tAI[codonName])
        totalCodonsProduct = float(np.prod(codonTaiTotalValues))
        totalCodonLength = float(len(totalCodons))
        totalTAI = totalCodonsProduct**(1/totalCodonLength)
        totalTAIs[shortID] = totalTAI

    # cumulativeAveragePositionalTAIs looks like this:
    # {0: 0.305623, 1: 0.2414229829717958, 2: 0.285033520896916, 3: 0.28437184018441025, 4: 0.33140008753305467, ....

    print ("cumulativeAveragePositionalTAIs",cumulativeAveragePositionalTAIs)
    return firstSegmentTAIs,secondSegmentTAIs,totalTAIs,cumulativeAveragePositionalTAIs
        
####################
# GLOBAL VARIABLES #
####################

stops = ['TAA','TGA','TAG']

TICs = ['TAC','AAC','TTC','TAT','AAG','GAG','AAT','GAT','ATC','GAC','GTG']

codonMap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
            "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
            "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
            "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
            "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
            "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
            "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
            "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

residueMap = {"F":("TTT","TTC"), "L":("TTA","TTG","CTT","CTC","CTA","CTG"), "S":("TCT","TCC","TCA","TCG","AGT","AGC"),
              "Y":("TAT","TAC"), "*":("TAA","TGA","TAG"), "C":("TGT","TGC"), "W":("TGG",), "P":("CCT","CCC","CCA","CCG"),
              "H":("CAT","CAC"), "Q":("CAA","CAG"), "R":("CGT","CGC","CGA","CGG","AGA","AGG"), "I":("ATT","ATC","ATA"),
              "M":("ATG",), "T":("ACT","ACC","ACA","ACG"), "N":("AAT","AAC"), "K":("AAA","AAG"), "V":("GTT","GTC","GTA","GTG"),
              "A":("GCT","GCC","GCA","GCG"), "D":("GAT","GAC"), "E":("GAA","GAG"), "G":("GGT","GGC","GGA","GGG")}

aaType = {}
aaValue = {}
aaType["R"] = "+"
aaValue["R"] = 'Blues'
aaType["H"] = "+"
aaValue["H"] = 'Blues'
aaType["K"] = "+"
aaValue["K"] = 'Blues'
aaType["D"] = "-"
aaValue["D"] = 'Greens'
aaType["E"] = "-"
aaValue["E"] = 'Greens'
aaType["S"] = "P"
aaValue["S"] = 'Purples'
aaType["T"] = "P"
aaValue["T"] = 'Purples'
aaType["N"] = "P"
aaValue["N"] = 'Purples'
aaType["Q"] =  "P"
aaValue["Q"] = 'Purples'
aaType["Y"] =  "P"
aaValue["Y"] = 'Purples'
aaType["C"] = "N"
aaValue["C"] = 'Oranges'
aaType["G"] =  "N"
aaValue["G"] = 'Oranges'
aaType["P"] =  "N"
aaValue["P"] = 'Oranges'
aaType["A"] =  "N"
aaValue["A"] = 'Oranges'
aaType["I"] =  "N"
aaValue["I"] = 'Oranges'
aaType["L"] =  "N"
aaValue["L"] = 'Oranges'
aaType["M"] =  "N"
aaValue["M"] = 'Oranges'
aaType["F"] =  "N"
aaValue["F"] = 'Oranges'
aaType["W"] =  "N"
aaValue["W"] = 'Oranges'
aaType["V"] =  "N"
aaValue["V"] = 'Oranges'

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

########
# MAIN #
########

usage = "Usage: python " + sys.argv[0] + " <fastafile> <PDCUB score file> <CDS window size in nt (multiple of 3)>"
# The window size is half the total sequence length. For example, a window size of 150 nt means a sequence length of 100 codons.

if len(sys.argv) != 4:
    print(usage)
    sys.exit()
    
fastaFile = sys.argv[1]
pdcubScores = sys.argv[2]
windowSizeNt = sys.argv[3]

binSizeNt = int(windowSizeNt)
binCountNt = 2 # number of bins to examine
ntSequenceLength = binCountNt*binSizeNt # this is total examined length of sequence, counting by nucleotides
codonSequenceLength = ntSequenceLength/3 # this is total examined length of sequence, counting by codons

pdcubScoresDict = readPdcubScores(pdcubScores)
valid_mRNAs = processFastaFile(fastaFile)

firstSegmentTAIs,secondSegmentTAIs,totalTAIs,cumulativeAveragePositionalTAIs = calculateTAIs(valid_mRNAs)

# color widget for coloring residue groups in plots
colorWidget = cm.get_cmap('viridis')

# plot tAI of second segment against tAI of first segment
plt.figure(dpi=600)

colors = ['midnightblue','lime']
colorMap = LinearSegmentedColormap.from_list('dah', colors, N=100)

taiSegOneVals = []
taiSegTwoVals = []
taiTotalVals = []
cols = []
IDs = []
for ID in firstSegmentTAIs.keys():
    IDs.append(ID)
for ID in IDs:
    taiSegOneVals.append(float(firstSegmentTAIs[ID]))
    taiSegTwoVals.append(float(secondSegmentTAIs[ID]))
    taiTotalVals.append(float(totalTAIs[ID]))
    cols.append(float(pdcubScoresDict[ID]))
minCol = min(cols)
maxCol = max(cols)
minTai = min(taiTotalVals)
maxTai = max(taiTotalVals)

scatterPoints = list(zip(taiSegOneVals,taiSegTwoVals,taiTotalVals,cols,IDs))
scatterPoints.sort(key=lambda x:x[3], reverse=False) # True means low PDCUB values will be drawn last and therefore be on top of plot
taiSegOneVals,taiSegTwoVals,taiTotalVals,cols,IDs = zip(*scatterPoints)

plt.scatter(taiSegOneVals,taiSegTwoVals,c=cols,s=5,cmap=colorMap,edgecolor='none',alpha=0.8)
slope, intercept, r_value, p_value, std_err = stats.linregress(taiSegOneVals,taiSegTwoVals)
r_squared = r_value*r_value
def bestfitGraph(regressionFormula, x_range):
    x = np.array(x_range)
    y = regressionFormula(x)
    # plt.plot(x,y,c='blue',label='$R^2$={:.2}'.format(r_squared),lw=1.2,linestyle='dashed')
    plt.plot(x,y,c='blue',label='{:.2}'.format(slope)+'x + '+'{:.2}'.format(intercept)+', $R^2$={:.2}'.format(r_squared),lw=1.2,linestyle='dashed')
def regressionFormula(x):
    return x*slope + intercept
bestfitGraph(regressionFormula, [x/10.0 for x in range(2,7,1)])
plt.colorbar(label='PDCUB score')
#for index,codonName in enumerate(RA2codons):
#    plt.annotate(codonName,(RA1values[index],RA2values[index]),fontsize=5)
#    if codonName in TICs:
#        plt.annotate(codonName,(RA1values[index],RA2values[index]),fontsize=5,color='red')
plt.plot([0.2,0.6],[0.2,0.6],label="y = x",c='black',lw=0.3,linestyle='dashed')
plt.xlim([0.2,0.6])
plt.ylim([0.2,0.6])
plt.legend(framealpha=0.5)
plt.xlabel("tAI (0 to " + str(binSizeNt) + " nt)")
plt.ylabel("tAI (" + str(binSizeNt) + " to " + str(binSizeNt*binCountNt) + " nt)")
plt.title("position-dependent tAI")
plt.savefig("posDepTAI_" + str(binSizeNt)  + "bin.pdf")
plt.savefig("posDepTAI_" + str(binSizeNt)  + "bin.png")
plt.close()

# plot tAI against PDCUB score
plt.figure(dpi=600)
plt.scatter(cols,taiTotalVals,c='black',s=5,edgecolor='none',alpha=0.8)
slope, intercept, r_value, p_value, std_err = stats.linregress(cols,taiTotalVals)
r_squared = r_value*r_value
bestfitGraph(regressionFormula, range(int(minCol),int(maxCol)))
plt.xlim([minCol,maxCol])
plt.ylim([minTai,maxTai])
plt.legend(framealpha=0.5)
plt.xlabel('PDCUB score')
plt.ylabel('local tAI of first ' + str(codonSequenceLength) + ' codons')
plt.title('Comparison of PDCUB and local tAI')
plt.savefig('taiVsPdcub.' + str(codonSequenceLength) + 'codons.pdf')
plt.savefig('taiVsPdcub.' + str(codonSequenceLength) + 'codons.png')
plt.close()
