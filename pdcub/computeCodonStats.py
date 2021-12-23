#Accepts a Gene FASTA file and calculates bin-based codon percentages, z-scores, p-values, and log likelihoods.
#Outputs each result in a TSV file.

#Author: Kaavya Subramanian




import sys, re
import math
from Bio import SeqIO
import statistics
import csv
import scipy.stats as st


stops = ['TAA','TGA','TAG']
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





#######################################################################################
#Class that represents a bin that the transcript is divided into. Holds information
#like the number of codons in that bin, and the zscore. Each codon has a list of bin
#objects that encompass the entire transcript.
#######################################################################################

class binObject:

    def __init__(self, nameofBin,codonCounts, zscore,percent):
            self.nameofBin = nameofBin      #useful when creating the heat map
            self.codonCounts = codonCounts
            self.zscore = zscore
            self.percent = percent

    def calcZscore(self, mean, sd):
        if sd == 0:     #when there are no codons in the transcript at that position
            self.zscore = 0
        else:
            self.zscore = (self.percent - mean) / sd


    def calcpercent(self,totalCounts):
        #print(totalCounts)
        self.percent = self.codonCounts/float(totalCounts)
        #try:
         #   self.percent = self.codonCounts/float(totalCounts)
          #  print(self, self.codonCounts, float(totalCounts))
        #except ZeroDivisionError:
         #   print(self.codonCounts, float(totalCounts))

#######################################################################
#Creates a list of Bin objects that every codon will use.
#######################################################################

def createBinList(binLength, maxLength):
    numberofBins = maxLength / binLength
    

    if maxLength % binLength != 0:
        numberofBins += 1   #in case it isn't a clean divide, then need to round up

    numberofBins = int(numberofBins)

    binList = []
    

    for x in range(0, numberofBins):
        
        binStartPos = (x*binLength) + 1
        binEndPos = (x+1) * binLength
        nameofBin = "%d - %d" % (binStartPos, binEndPos)
        newBin = binObject(nameofBin,0,0,0) #initialize zscore and codon count values as 0
        #print(newBin.codonCounts)
        binList.append(newBin)
        #print(binList)

    return binList



####################################################################
#Creates a dictionary where the key is the codon and is associated
#with a list of bins.
###################################################################

def createCodonDictionary(binLength, maxLength):
    codonVals = {}
    for codon in codonMap.keys():
        codonVals[codon] = createBinList(binLength,maxLength)

    return codonVals


#Calculates the percents for all codons in the dictionary.

def calcAllPercents(codonVals,allBinValues):
    for key in codonVals.keys():
        #print(codonVals)
        n = 0
        for x in codonVals[key]:
            #print(x)
            x.calcpercent(allBinValues[n])
            n += 1
#Calculates the z-score for all codons in the dictionary

def calczscores(codonVals):
    for key in codonVals.keys():
        newData = []      #each codon has a unique mean and sd.
        for x in codonVals[key]:
            newData.append(x.percent)

        mean,sd = 0.0,0.0
        if len(newData) > 1:    #there must be at least two values to compute the mean and sd.
            mean = statistics.mean(newData)
            sd = statistics.stdev(newData)
        for x in codonVals[key]:
            x.calcZscore(mean,sd)

#Sequence Classes and Functions. Returns the CDS of a given transcript based on headerline information.

def parseDefline(record):
    #print(record.seq)
    #match  = re.search('CDS:(\d+)-(\d+)', record.id)
    #start = int(match.group(1))
    #end = int(match.group(2))
    #rnaSeq = record.seq.transcribe()
    #cds = rnaSeq[start-1:end]
    cds = record.seq
    #cdsRange = len(cds)
    return cds



#Checks to see if a CDS has any anamolies. Only valid CDS transcripts are returned
def checkGivenCodingSequence(geneTranscript):
    startCodon = "ATG"

    transcriptStartCodon = geneTranscript[0:3]
    transcriptStopCodon = geneTranscript[len(geneTranscript)-3:len(geneTranscript)]


    if (transcriptStartCodon == startCodon) and (transcriptStopCodon in stops) and (len(geneTranscript) % 3 == 0):
        return True
    else:
        return False

#Parses transcripts in FASTA file and extracts CDS. Filters out valid CDSs.
def returnValidSequences(fastaFile, allCodonValues):

    print("Extracting valid sequences from {}.".format(fastaFile))
    longestCDS = 0
    validSequences = []

    for record in SeqIO.parse(fastaFile, "fasta"):
        currentSequence = parseDefline(record)
        #print(currentSequence)
        check = checkGivenCodingSequence(currentSequence)
        #print(check)
        if check:

            #10/30/2020: get codon global counts:
            rnaSeq = record.seq.transcribe()
            for key in codonMap.keys():
                allCodonValues[key] += rnaSeq.count(key)

            currentSequence = currentSequence[3:]  #exclude start codon
            validSequences.append(currentSequence)
            if len(currentSequence) >= longestCDS:
                longestCDS = len(currentSequence)  #Keeps track of the longest CDS to create the necessary amount of bins.
            else:
                 continue

    print("Extracted {} sequences.".format(len(validSequences)))
    return validSequences, longestCDS

binValues = int(sys.argv[3]) 

def countCodonsintoBins(validSequences,maxLength,codonVals,allBinValues,allCodonValues):
    num = 0     #Keeps track of the nth transcript read.
    #print(validSequences, len(validSequences))
    for x in validSequences:
        #print(x)
        num += 1
        sys.stdout.write("Transcript %d/%d    \r" %(num,len(validSequences)))
        sys.stdout.flush()
        for y in range(0,len(x),3):
            #print(y)
            if y < maxLength:
                codon = x[y:y+3]
                binIndex = y//binValues
                codonVals[codon][binIndex].codonCounts += 1
                allCodonValues[codon] += 1
                allBinValues[binIndex] += 1
                
                #print(binIndex, allBinValues[binIndex], codon, codonVals[codon][binIndex].codonCounts, allCodonValues[codon], allBinValues[binIndex])

def createTSV(filename,data,getBins):
   # with open(filename,'wb') as f_output:
    with open(filename,'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter = '\t')
        tsv_output.writerow(getBins)
        for row in data:
            tsv_output.writerow(row)

def main():

    usage = "Usage: " + sys.argv[0] + " <Gencode FASTA>" + "<results Directory>" + "<bin size (nt)>" + "<CDS cutoff length (nt)>"
    if len(sys.argv) != 5:
        print(usage)
        sys.exit()

    fastaFile = sys.argv[1]
    resultsDir = sys.argv[2]
    binValues = int(sys.argv[3])   #predetermined value
    maxLength = int(sys.argv[4])

    print("Starting.")

    #stops = ['TAA','TGA','TAG']
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




    print("Parsing FASTA file.")

    allCodonValues = {}             #Dictionary that keeps track of the total number of codons in all transcripts.
    for key in codonMap.keys():
        allCodonValues[key] = 0

    validSequences, longestCdsLengths  =  returnValidSequences(fastaFile, allCodonValues)
    #Note on Validsequences: All sequences are modified to exclude start codon



    #Create Data Structures to hold codon values.
    codonVals = createCodonDictionary(binValues,maxLength)


    allBinValues = {}
    for b in range(len(codonVals['TTT'])):
        allBinValues[b] = 0



    print("Counting codons into their bins.")
    #print(countCodonsintoBins(validSequences,maxLength,codonVals,allBinValues,allCodonValues))
    countCodonsintoBins(validSequences,maxLength,codonVals,allBinValues,allCodonValues)



    print("Calculating percentages.")
    calcAllPercents(codonVals,allBinValues)

    print("Calculating z-score values.")
    calczscores(codonVals)

    #Creates list of Codon names and bin ranges as labels.
    codonNames = []
    for x in residueMap.keys():
        if x != "*":
            for y in residueMap[x]:
                codonNames.append(y)

    getBins = []
    getBins.append('Codon')
    for x in codonVals['TTT']:
        getBins.append(x.nameofBin)


        #Creates arrays that can be read into TSVs.
    zScoreData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in codonVals[x]:
            val = 0 + y.zscore
            newList.append(val)
            zScoreData.append(newList)


    print("Calculating z-squared values.")

    zSquareData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in codonVals[x]:
            val = y.zscore ** 2
            newList.append(val)
            zSquareData.append(newList)


    print("Calculating p-values.")
    #Uses Zscore data to do it since it is essentially a single function change.
    pvalData = []
    for codonRow in zScoreData:
        newList = []
        newList.append(codonRow[0])
        z_counter = 1
        for z_counter in range(1,len(codonRow)):
            p = st.norm.sf(codonRow[z_counter])    #Finds upper tail probability.
            newList.append(p)
        pvalData.append(newList)

    print("Counting total numbers of each codon per bin.")
    codonCountData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in codonVals[x]:
            val = 0 + y.codonCounts
            newList.append(val)
        codonCountData.append(newList)

    totalCodons = 0
    for x in allCodonValues:
        totalCodons += allCodonValues[x]



    print("Calculating log-likelihood score.")
    logLikelihoodData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in codonVals[x]:
            z = allCodonValues[x]/float(totalCodons)
            weight = float(y.percent)/(z)
            #val = math.log(weight,2)
            try:
                val = math.log(weight,2)
                #print("pass", weight)
            except ValueError:
                #print("val error", weight)
                val = 0
            newList.append(val)
        logLikelihoodData.append(newList)


    percentData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in codonVals[x]:
            val = 0 + y.percent
            newList.append(val)
        percentData.append(newList)


    print("Saving results.")
    createTSV('{}/codonCounts.tsv'.format(resultsDir),codonCountData,getBins)
    createTSV('{}/zscores.tsv'.format(resultsDir),zScoreData,getBins)
    createTSV('{}/pvalues.tsv'.format(resultsDir),pvalData,getBins)
    createTSV('{}/zsquaredvals.tsv'.format(resultsDir), zSquareData,getBins)
    createTSV('{}/percentdata.tsv'.format(resultsDir), percentData, getBins)
    createTSV('{}/logLikelihoodData_{}.tsv'.format(resultsDir,binValues),logLikelihoodData, getBins)

if __name__ == '__main__':
    main()
