#This is the main Thesis file. It does't work from top to down with one run.
#I usually run per piece to get what I need.
#In this file I am parsing at first the 2 GFF3 files i have (Nagalakshmi and Yassour), that give information
#about the sequence (eg 5UTR-start position, end Position, Chromosome etc)
#Later On I constract different Features and I do the Regression

from __future__ import with_statement

import copy
import csv
import gzip
import math
import matplotlib.pyplot as plt
import pdb
import plotly.plotly as py
import numpy as np
np.set_printoptions(threshold='nan')
import scipy
import time
import urllib

from collections import namedtuple, Counter
from copy import deepcopy
from plotly.graph_objs import *
py.sign_in('Eftychia', '2puhmq6aj8')
from pylab import *
from scipy.stats.stats import pearsonr
from sklearn import datasets, linear_model, cross_validation, metrics, clone, gaussian_process, svm, preprocessing
from sklearn.cross_validation import KFold, cross_val_score, StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, ExtraTreesClassifier, AdaBoostClassifier
from sklearn.externals.six.moves import xrange
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import Pipeline

#Import My Own Functions
from myfunctions import parseGFFAttributes, parseGFF3, find_kmers, romanToNumeric, shortenSequence, findDuplcicates, searchIfInColumn, plotWith4Subplots, plotTwoScales, findBasesFrequency

#filename = "Nagalakshmi_2008_UTRs.gff3"
#filename = "Yassour_2009_UTRs.gff3"

#######################################################################
#######Chromosomes
'''
chromosomes = []
chromosomes.append('') # i want to have empty string in 0 position.

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip() # rstrip trims the white spaces from the right
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq)) #yield is like return, except the function will return a generator

for i in range(1,18):
    with open("chr"+`i`+".fsa") as fp:
        for name, seq in read_fasta(fp):
            #print(name, seq)
            chromosomes.append(seq)

#Now I have an array of Chromosomes that has all the 17 chromosomes 0,1....17
#print chromosomes[16]
#print "--------------------------------------------------------------------------------------------------end of chromosomes"
'''
#######################################################################
'''
#Open also the yeast names file withcsvFile = 'newDataset.csv' the translation initiation rates
# csvFile = 'newDataset.csv'
csvFile = 'ciandrini.csv'
genes = [] #has the name of the genes
InitiationRates = [] # has the translation initiation rates-

with open(csvFile,'rb') as csvfile: #read binary file.
    data = list(csv.reader(csvfile, delimiter=','))
    #print data
    teller = 0
    for row in data:
        teller += 1
        if len(row) > 1:
            #print ', '.join(row)
            if "-" in row[0]:
                row[0] = row[0][:-2]
            if "#" in row[0] or "#" in row[1]:
                #print row
                continue
            genes.append(row[0])
            InitiationRates.append(math.log(float(row[1]))) #base e
        else:
            print "row; ",teller,row

# np.save("InitRatesCiandrini",InitiationRates)
# exit()

# print genes
# print InitiationRates
# print "--------------------------------------------------------------------------------------------------end of csv file."
'''
########################################################################
#Parse gff3 file.
'''


#Initialized GeneInfo named tuple. Note: namedtuple is immutable(ametavlhtos)
gffInfoFields = ["source", "filename", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)



#end of functions

#------------------------------------------------------------------------------------------

if __name__ == "__main__":

    existingGenes = []
    actualArray   = []
    sequence = []
    upstreamsequences = []

    recordCount = 0
    for record in parseGFF3("Nagalakshmi_2008_UTRs.gff3", gffInfoFields, GFFRecord):

        #print record

        #Access attributes like this: my_strand = record.strand
        # print "my_strands_name is ", record.attributes['Name'][:-5]

        #Roman to Numerical to find which Chromosome
        recSource = romanToNumeric(record.source)
        
        name = record.attributes['Name'][:-5]
        wholeName = record.attributes['Name']
        # print "wholeName", wholeName
        
        #Some genes have name like this YPR036W-A but in order to check if they are in the csv file i need only the first part.
        pos = record.attributes['Name'][:-5].find('-')

        if pos != -1 :
            foundName = record.attributes['Name'][:pos]
        else:
            foundName = name #record.attributes['Name'][:-5]
        
        #Now check if this name exists also in the csv file.
        if (foundName in genes) and ('5UTR' in wholeName):

            fiveUTRStart = shortenSequence(record.start,record.end) # I want only the last 100bp upsream
            fiveUTREnd = record.end + 40 # I want to get the First 40bp of the Coding Sequence
            sequence = chromosomes[recSource][fiveUTRStart:fiveUTREnd+1]
            #print "sequence existing", sequence
            sequenceLength = len(sequence)
            #print "sequence %s and %d" % (sequence,sequenceLength)

            #For the MFE and Length I want the whole 5UTR
            #print "-->>", record.start, record.end
            seq = chromosomes[recSource][record.start:record.end+1]
            sequenceLe = len(seq)
            #print sequenceLe

            dimers,trimers,dimercounts,trimercounts =find_kmers(sequence)

            #Bases Frequencies
            Afrequency, Tfrequency, Gfrequency, Cfrequency = findBasesFrequency(sequence)
            # print Afrequency, Tfrequency, Gfrequency, Cfrequency
            
            #baseCounter = Counter(sequence)
            #Prints : baseCounter Counter({'A': 64, 'C': 27, 'G': 26, 'T': 23})


            #append them in an array
            geneindex = genes.index(foundName)
            
            #if (InitiationRates[geneindex] < -9) and (InitiationRates[geneindex] > -12) : # For Alexey's Data I m interested in [-12,-9]


            upstreamsequences.append([seq, sequenceLe])
            #print "upstreamsequences", upstreamsequences
            

            # print "InitiationRates[geneindex]", InitiationRates[geneindex]
            # print "YES", geneindex

            existingGenes.append([foundName,recSource,sequenceLength,Afrequency,Tfrequency,Gfrequency,Cfrequency,InitiationRates[geneindex]]) #this array has all the matching genes 

            instanceArray = [sequence,recSource,sequenceLength,Afrequency,Tfrequency,Gfrequency,Cfrequency]

            maxNoDimers = float(sequenceLength-1)
            maxNoTrimers = float(sequenceLength-2)
            for dim in dimercounts:
                if maxNoDimers != 0:
                    instanceArray.append(dim/maxNoDimers)
                else:
                    instanceArray.append(dim)
            for trim in trimercounts:
                if maxNoTrimers != 0:
                    instanceArray.append(trim/maxNoTrimers)
                else:
                    instanceArray.append(trim)

            instanceArray.append(InitiationRates[geneindex])
            #print instanceArray
            actualArray.append(instanceArray) #this has the thresholded 

        recordCount += 1


    # print "The Nagalakshmi file has %d rows" % recordCount
    # print "The existingGenes compared with Nagalakshmi has %d rows" % len(existingGenes)
    # print "the actualArray has %d rows" % len(actualArray)
    # print "existingGenes", existingGenes
    # print "actualArray", actualArray
    transExistingGenes = np.array(existingGenes).T[0]
#    print np.array(existingGenes).T[0]
    
########and again for the second file.################################################## 2nd File
    recordCount = 0
    for record in parseGFF3("Yassour_2009_UTRs.gff3", gffInfoFields, GFFRecord):

        #Roman to Numerical to find which Chromosome
        recSource = romanToNumeric(record.source)

        name = record.attributes['Name'][:-5]
        wholeName = record.attributes['Name']
        #Some genes have name like this YPR036W-A but in order to check if they are in the csv file i need only the first part.
        pos = record.attributes['Name'][:-5].find('-')

        if pos != -1 :
            foundName = record.attributes['Name'][:pos]
        else: # if -1 means there is not
            foundName = name # record.attributes['Name'][:-5]

        #Now check if this name exists also in the eftychia.csv file.
        if (foundName in genes) and (foundName not in transExistingGenes) and ('5UTR' in wholeName):

            #COunt the number of A, T, C, G
            fiveUTRStart = shortenSequence(record.start,record.end) # I want only the last 100bp upsream
            fiveUTREnd = record.end + 40 # I want to get the First 40bp of the Coding Sequence
            sequence = chromosomes[recSource][fiveUTRStart:fiveUTREnd+1]
            #print "sequence", sequence
            sequenceLength = len(sequence)

            #For the MFE and Length I want the whole 5UTR
            seq = chromosomes[recSource][record.start:record.end+1]
            sequenceLe = len(seq)

            dimers,trimers,dimercounts,trimercounts =find_kmers(sequence)


            #Bases Frequencies
            Afrequency, Tfrequency, Gfrequency, Cfrequency = findBasesFrequency(sequence)
            # print Afrequency, Tfrequency, Gfrequency, Cfrequency
            
            #baseCounter = Counter(sequence)

            geneindex = genes.index(foundName)
            # if (InitiationRates[geneindex] < -9) and (InitiationRates[geneindex] > -12) :
                

            upstreamsequences.append([seq, sequenceLe])

            existingGenes.append([foundName,recSource,sequenceLength,Afrequency,Tfrequency,Gfrequency,Cfrequency,InitiationRates[geneindex]]) #whole array

            instanceArray = [sequence,recSource,sequenceLength,Afrequency,Tfrequency,Gfrequency,Cfrequency]

            maxNoDimers = float(sequenceLength-1)
            maxNoTrimers = float(sequenceLength-2)
            for dim in dimercounts:
                if maxNoDimers != 0:
                    instanceArray.append(dim/maxNoDimers)
                else:
                    instanceArray.append(dim)
            for trim in trimercounts:
                if maxNoTrimers != 0:
                    instanceArray.append(trim/maxNoTrimers)
                else:
                    instanceArray.append(trim)

            instanceArray.append(InitiationRates[geneindex])
            #print instanceArray
            actualArray.append(instanceArray) #this has the genes

            recordCount += 1
    # print "The Yassour file has %d rows" % recordCount

    # print "The existingGenes has now %d rows" % len(existingGenes)
#    print existingGenes
    print "the actualArray has now %d rows" % len(actualArray)

# np.save("existingGenesAlexey", existingGenes)
# np.save("actualArrayAlexey", actualArray)
np.save("Ciandrini/existingGenes", existingGenes)
np.save("Ciandrini/actualArray", actualArray)

upstreamsequences = np.array(upstreamsequences)
print "upstreamsequences", upstreamsequences.shape
np.save('Ciandrini/upstreamsequences', upstreamsequences)
'''
upstreamsequences = np.load('Ciandrini/upstreamsequences.npy')
existingGenes     = np.load('Ciandrini/existingGenes.npy')
actualArray       = np.load('Ciandrini/actualArray.npy')
print "upstreamsequences loades", upstreamsequences.shape
print "existingGenes loaded", existingGenes.shape
print "actualArray loaded", actualArray.shape

print ""

#print "--------------------------------------------------------------------------------------------------end of parsing gff3 files"

#######################################################################################################################################################
#######################################################################################################################################################
# #This is only for Alexey's data!!
# extraAlexey = np.load('Gritsenko/extraAlexey.npy') #For the extra Numbers
# print "extraAlexey", extraAlexey.shape, "this is array with mRNA count, RibosomeCount, Fitness, Segments"
# #print "extraAlexey", extraAlexey[1]

# actualArray   = np.array(actualArray)    # [Sequence, Chromosome, Length, A, T, G, C, 2mers, 3mers, Initiation Rates]
# print "actualArray", actualArray.shape      # actualArray (3409, 88) and it Contains Sequences

# existingGenes = np.array(existingGenes)  # [foundName,Chromosome,Length,Afrequency,Tfrequency,Gfrequency,Cfrequency,InitiationRates]
# # print "existingGenes", existingGenes[0]
# print "existingGenes", existingGenes.shape  # existingGenes (3409, 8)

##STRESS##
#Load the Genes that are related to Stress
stressedGenes = np.load('Ciandrini/stressedGenes.npy')
print "stressedGenes loaded", stressedGenes.shape
#stressedGenes ['YAL062W' '0']

actualArrayWithNames = np.vstack((existingGenes[:,0],actualArray.T)).T    #(3409, 89) [Name, Sequence, Chromosome,........., InitRates]
# np.save('Ciandrini/actualArrayWithNames', actualArrayWithNames)

actualArrayWithNames = np.load('Ciandrini/actualArrayWithNames.npy')
print "actualArrayWithNames loaded", actualArrayWithNames.shape
print ""


#Add One more column of 0s
instanceArray = actualArrayWithNames[:,0:-1]
noRows = len(actualArrayWithNames)
print"noRows", noRows
z = np.zeros(noRows)
instanceArray = np.vstack((instanceArray.T,z)).T

#Now Change tHe Length into the Length of the 5UTR && Add The Stress Feature
for i in range(noRows):
    
    #They are in the same order
    instanceArray[i][3] = upstreamsequences[i][1]

    #Now add Stress, they are also in the same order
    index = np.where(instanceArray[i][0]==stressedGenes[:,0])
    index = index[0][0]
    # print index

    instanceArray[i][-1] = stressedGenes[index][1]
print "Stress added as feature instanceArray:", instanceArray.shape

mfe = np.load('Ciandrini/mfe.npy')
print "Mfe Loaded", mfe.shape

instanceArray = np.vstack((instanceArray.T, mfe)).T

"Mfe added as feature", instanceArray.shape

instanceArray = np.vstack((instanceArray.T,actualArrayWithNames[:,-1])).T #InitRates
featureArray = instanceArray #(3409, 90)
print "featureArray", featureArray.shape
np.save('Ciandrini/featureArray', featureArray)
exit()




#Correlation
pearsonsCorrelations = []
spearmanCorrelations = []
# rawFeatures = actualArrayWithNames[:,2:-1]
# rawTarget = actualArrayWithNames[:,-1]
Features = actualArrayWithNames[:,2:-1].astype(float)
Target = actualArrayWithNames[:,-1].astype(float)
print "Features", Features.shape

############################################################################################################################################
##########################################################################################################################################
'''
sortedArray = []
initRates= actualArrayWithNames[:,-1].astype(float)
sortedArray = actualArrayWithNames[initRates.argsort()] # Sort the array by InitiationRate
sortedArray = sortedArray[:,1:]
'''

print ""
print "actualArrayWithNames", actualArrayWithNames.shape
print ""
actualArray = actualArrayWithNames[:,1:] #Update the Values of the ActualArray!!!!!!!!!!!!!
# np.save('Gritsenko/forCE', actualArrayWithNames)
# exit()

#################################################################################################################################################
################## WRITE a file, save an array/list into a file.
'''
# ~Save the sequences in a fasta file to use it for the energy.
mySecs = actualArray[:,0]
#f = open('mySecsCiandrini.fa','w')
f = open('Stress/mySecsAlexey.fa','w')
for item in mySecs:
    f.write('>sec\n%s\n' % item)
f.close
exit()
'''
#~~~Add the extra feature to the main array~~~
print "actualArray", actualArray.shape #,actualArray[0]

print "instanceArray = deepcopy(actualArray[:,1:-1])"
print ""
instanceArray = []
instanceArray = deepcopy(actualArray[:,0:-1]) # Copy all features 

sumEsetSequences = np.load("Gritsenko/sumEsetSequencesAlexey.npy")
print "Summation of the Conditional Entropy values is done: sumEsetSequences=", sumEsetSequences.shape

instanceArray = np.vstack( (instanceArray.T, sumEsetSequences ) )    # ADD the SUmmation of Entropies Feature
print "Summations Added", instanceArray.shape
# instanceArray[0,:] those are all the lengths
print "this is the first row with all the features including the SuMmAtIoN", instanceArray[:,0]
print ""


mfe = [] #Minimum Free energy
# f = open('mfesciandrini.txt','r')
f = open('Stress/mfesalexey.txt','r')
for line in f:
    mfe.append(line)
mfe = np.array(mfe)
mfe = mfe.astype('float')
print "mfe", mfe.shape


instanceArray = np.vstack( (instanceArray, mfe) ) # ADD the MFE as feature
print ""
print "~~~ MFE added in the Feature array: instanceArray", instanceArray.shape
# print "instanceArray[1]", instanceArray[1], "Those are the lengths"
# print ""
# print "instanceArray[:,0]", instanceArray[:,0] , "This is the first row with all the features including mfe"
print ""

#Pre-work before adding them to the array
extraNumbers = np.zeros((len(actualArrayWithNames),4)) #Create 0 containing columns
for j in range(len(extraAlexey)):
    for i in range(len(actualArrayWithNames)):
        if actualArrayWithNames[i][0] == extraAlexey[j][0]:
            #print actualArrayWithNames[i][0], extraAlexey[j]
            extraNumbers[i][0] = extraAlexey[j][2]   # mRNA read count
            extraNumbers[i][1] = extraAlexey[j][3]   # Ribosome read count
            extraNumbers[i][2] = extraAlexey[j][4]   # Fitness
            extraNumbers[i][3] = extraAlexey[j][5]   # No of Segments

print "extraNumbers", extraNumbers.shape, extraNumbers[0:5]
print ""
print "-->instanceArray", instanceArray.shape, extraNumbers.shape
featureArray = np.vstack((instanceArray, extraNumbers[:,0])) #mRNA counts
featureArray = np.vstack((featureArray, extraNumbers[:,1])) #Ribosome Counts
featureArray = np.vstack((featureArray, extraNumbers[:,2])) #Fitness
featureArray = np.vstack((featureArray, extraNumbers[:,3])) #No of Segments

#Add AvgFitness as a feature
b = np.zeros((len(extraNumbers)))
featureArray = np.vstack((featureArray, b))
featureArray = featureArray.T
for i in range(len(featureArray)):
    featureArray[i][-1] = (featureArray[i][-3].astype(float)/featureArray[i][-2].astype(float))
    # print featureArray[i][-1]
featureArray = featureArray.T



print "Added the Extra Numbers from Alexey"
print "featureArray", featureArray.shape
print ""


featureArray = np.vstack((featureArray, actualArray[:,-1])) #add the initiation rates

featureArray = featureArray.T #.astype(float) #this is the correct shape
print "Initiation Rates added: featureArray", featureArray.shape
#print "featureArray", featureArray[0:5, :]
print "--"

temp = np.array(featureArray[:,-2],dtype=float)
where_are_NaNs = isnan(temp) #Remove the nans!
temp[where_are_NaNs] = 0
featureArray[:,-2] = temp

np.save("Gritsenko/featureArray", featureArray)
exit()
# sortedar = []
# initRates= featureArray[:,-1].astype(float)
# sortedar = featureArray[initRates.argsort()] # Sort the array by InitiationRate
# # sortedar = sortedar[:,1:]
# #Plot Summations
# y1 = sortedar[:,88]
# y2 = sortedar[:,-1]
# x  = np.arange(len(y1))
# plotTwoScales(x, 'samples', y1, 'Summations', y2, 'InitRates', 'SummedValues from CE', 'Initiation Rates (log)', 'Gritsenko/plotSummations.png')
# exit()


# Now all The Features Are Collected.
###################################################################################################################################################
'''
# Pre Selection OF features
# Depending on their Correlation

X = featureArray[:,0:-1]
y= featureArray[:,-1]

#Calculate the Correlation at each fold
pearsonsCorrelations = []
spearmanCorrelations = []

for i in range(len(X.T)): #per column.   # I should use Features and not rawFeatures
    pC = scipy.stats.pearsonr(X[:,i], y)
    pearsonsCorrelations.append(pC)
    # prawrint pC
    sC = scipy.stats.spearmanr(X[:,i], y)
    spearmanCorrelations.append(sC)
    # print sC

pearsonsCorrelations = np.array(pearsonsCorrelations)
spearmanCorrelations = np.array(spearmanCorrelations)
# print "pearsonsCorrelations", pearsonsCorrelations.shape
# print "spearmanCorrelations", spearmanCorrelations.shape

h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
z = spearmanCorrelations[:,0] # Spearman

#Use Spearman Correlation
selectedFeats = []
for i in range(len(z)):
    if math.fabs(z[i]) > 0.04:
        print "math.fabs(z[i])", i, z[i],  math.fabs(z[i])
        selectedFeats.append(i)
print "selectedFeats", selectedFeats

'''
print "featureArray before regression", featureArray.shape
# print featureArray[0]
# print "Print Per column the featureArray"
# for i in range(len(featureArray[0].T)):
#     print i, featureArray[0][i]
# print "-->", featureArray[0][93]


#Do Gene Selection
geneSelection = []
count = 0
# elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],80) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],90):
#         featureArray[seq,-1] = 9
#     elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],90):
#         featureArray[seq,-1] = 10

#Sort the array according to avg Fitness
avgFitness=featureArray[:,-2].astype(float)
sortdFeatsArry = featureArray[avgFitness.argsort()]
#The sorting goes from Small to Big, I am interested in the 10% Highest which is last!

print "np.percentile(featureArray[:,-2],90", np.percentile(sortdFeatsArry[:,-2],90)
threshold = np.percentile(sortdFeatsArry[:,-2],90)
for i in range(len(featureArray)):
    if featureArray[i][-2] > threshold:
        geneSelection.append(featureArray[i])
    else:
        count +=1


# for i in range(len(featureArray)):
#     if featureArray[i][-2] < 0:
#         #print i, featureArray[i][-2]
#         count += 1
#     else:
#         geneSelection.append(featureArray[i])
geneSelection = np.array(geneSelection)
print "geneSelection", geneSelection.shape
print "not selected genes", count

X = geneSelection[:,1:-6] #featureArray[:,selectedFeats] # i dont want the Chromosome Number as a feature which is the 1st
y= geneSelection[:,-1]
print "X", X.shape
print "y", y.shape


skf = cross_validation.KFold(len(y),n_folds=5)
# skf = cross_validation.StratifiedKFold(b,n_folds=5)

#Do Random Forest
print "RandomForest Start"
print ""
allScores = []
fold = 1
for train_index, test_index in skf:
    # print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    rndmForest = RandomForestRegressor(n_estimators=4000)
    print "X_train",X_train.shape
    print "y_train",y_train.shape

    print "Start training Random Forest"
    start_time_training = time.time()
    rndmForest.fit(X_train, y_train)
    print "Done training in ",time.time() - start_time_training, "seconds"
    
    score = rndmForest.score(X_test, y_test)
    allScores.append(score)
    y_predictTest = rndmForest.predict(X_test)
    y_predictTrain = rndmForest.predict(X_train)

    ys = np.vstack([y_train , y_predictTrain]).T
    ys = ys[ys[:, 0].argsort()]
    
    y_train = ys[:,0]
    y_predictTrain = ys[:,1]

    ys = np.vstack([y_test , y_predictTest]).T
    ys = ys[ys[:, 0].argsort()]
    
    y_test = ys[:,0]
    y_predictTest = ys[:,1]

    #Images
    plt.close('all')
    plt.figure()
    plt.subplots_adjust(hspace=.5, wspace=.5)
    plt.subplot(1,2,1)
    plt.title('rndmForest train')
    allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
    plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=3, label='train')
    plt.hold('on')
    plt.scatter(allSampleIndexes, y_predictTrain, c='r', edgecolor='none', s=3, label='prediction')
    plt.xlabel('#of rows')
    plt.ylabel('initiation rates (log)')
    art = []
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)

    plt.subplot(1,2,2)
    allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
    plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=3, label='test')
    plt.hold('on')
    plt.scatter(allSampleIndexes, y_predictTest, c='r', edgecolor='none', s=3, label='prediction')
    plt.xlabel('#of rows')
    plt.ylabel('initiation rates (log)')
    plt.title('rndmForest test')
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)
    plt.suptitle('rndmForest Gene Selection'+`fold`+'Fold')
    plt.savefig('GeneSelection/rndmForest'+ `csvFile` +'Fold '+`fold`+'.png', additional_artists=art, bbox_inches='tight')
    # plt.show()
    plt.clf()

    #Calculate the Correlation at each fold
    pearsonsCorrelations = []
    spearmanCorrelations = []
    
    for i in range(len(X_train.T)): #per column.   # I should use Features and not rawFeatures
        pC = scipy.stats.pearsonr(X_train[:,i], y_train)
        pearsonsCorrelations.append(pC)
        # prawrint pC
        sC = scipy.stats.spearmanr(X_train[:,i], y_train)
        spearmanCorrelations.append(sC)
        # print sC

    pearsonsCorrelations = np.array(pearsonsCorrelations)
    spearmanCorrelations = np.array(spearmanCorrelations)
    # print "pearsonsCorrelations", pearsonsCorrelations.shape
    # print "spearmanCorrelations", spearmanCorrelations.shape

    h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
    z = spearmanCorrelations[:,0] # Spearman
    x = np.arange(len(X_train.T)) # Features
    p = pearsonsCorrelations[:,1] # Pearson P values
    q = spearmanCorrelations[:,1] # Spearman P values

    plt.clf()
    ind = np.arange(len(X_train.T)) #width of a bar
    f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
    f.tight_layout()
    ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
    ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
    ax1.set_xticks(x)
    ax1.grid()
    ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
    ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
    ax2.set_title('P value')
    ax2.set_xticks(x)
    ax2.grid()
    f.suptitle('Gene Correlation at fold '+ `fold`)
    f.savefig('GeneSelection/RndmForestCorrelationFold'+ `fold`+'.svg')
    plt.clf()
    fold += 1

# print "Prediction of X is ", Y_predicted
allScores = np.array(allScores)
print "Random Forest Scores:", allScores
print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_predictTrain, y_train)))
print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_predictTest, y_test)))
exit()

###################################################################################################################################################
'''
#transform every sequence to an array of numbers
allNumberedSequences = []
for i in range(len(actualArray)): #pass from each sequence
    numberedSequence = []
    for j in range(maxLength): #pass from each position
        lenseq = int(actualArray[i][1]) # 1:Length, 0:Sequence
        if j < lenseq:
            if actualArray[i][0][lenseq-(j+1)] == 'A': #walk from right to left
                numberedSequence.append(1)
            elif actualArray[i][0][lenseq-(j+1)] == 'T':
                numberedSequence.append(2)
            elif actualArray[i][0][lenseq-(j+1)] == 'G':
                numberedSequence.append(3)
            elif actualArray[i][0][lenseq-(j+1)] == 'C':
                numberedSequence.append(4)
        else :
            numberedSequence.append(0)
    allNumberedSequences.append(numberedSequence)
np.save("allNumberedSequences", allNumberedSequences)
'''
allNumberedSequences = np.load("allNumberedSequences.npy")
print "Sequences are turned into numbers 1,2,3,4: allNumberedSequences=", allNumberedSequences.shape

################################ HEATMAPS - PLOTS ##########################################################################
#~~~~~Sorting~~~~~~~~~~~~
sortedFeaturesArray = []
print "--> Sorting, featureArray ", featureArray.shape
initRates=featureArray[:,-1].astype(float)
sortedFeaturesArray = featureArray[initRates.argsort()] # Sort the array by InitiationRate
print "sortedFeaturesArray:", sortedFeaturesArray.shape


temp = np.array(featureArray[:,-2],dtype=float)
maxValue=temp.max(axis=0) #axis=0:column: Indicates the no of Points/rows.
minValue=temp.min(axis=0)
print "maxValue of avg Fitness", maxValue
print "minValue of avg Fitness", minValue
exit()

# #Plot avg Fitness againt Init Rates
# x1 = sortedFeaturesArray[:,-2] # Avg. Fitness
# y  = sortedFeaturesArray[:,-1] # Initiation Rates
# x  = np.arange(len(y))
# plt.figure()
# plt.scatter(x, x1, c='r',s=1.5,lw = 0, label = 'Avg Fitness')
# plt.scatter(x, y, c='g',s=1,lw = 0, label = 'IntRates')
# plt.grid(True)
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=2, mode="expand", borderaxespad=0.)
# plt.savefig('avgFitnessPlot.png')


'''
plt.clf()
#~~Make a plot of Initiation Rates~~
plt.close('all')
plt.figure()
# plt.subplot(1,2,1)
y = sortedFeaturesArray[:,-1] #Initiation rates
x = np.arange(len(y))
plt.plot(x,y,'g', label='Alexey')
plt.xlabel("Indices")
plt.ylabel("Init Rates")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.grid(True)
plt.show()
exit()
plt.hold('on')
plt.subplot(1,2,2)
w = np.load("InitRatesCiandrini.npy")

x = np.arange(len(w))
plt.plot(x,w,'r', label='Ciandrini')
plt.xlabel("Indices")
plt.ylabel("Init Rates")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.grid(True)
# plt.show()
# exit()
plt.savefig("InitRatesBothAlexeyCiandrini.png")
plt.clf()
exit()

#~~Make a plot of Length~~
y = sortedFeaturesArray[:,0] #Length
z = sortedFeaturesArray[:,-1] #Initiation Rates
x = np.arange(len(y))
plt.figure()
plt.plot(x,y,'g', label='length')
plt.plot(x,z,'r', label='InitiationRates')
plt.xlabel("Indices")
plt.ylabel("Length")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
# plt.show()
plt.savefig("Plot of the Length for " + `csvFile` + ".png")
plt.clf()

#~~Make a plot of Summations~~
y = sortedFeaturesArray[:,-3] #sumEsetSequences
z = sortedFeaturesArray[:,-1] #Initiation Rates
print actualArray.shape
print y.shape
x = np.arange(len(y))
plt.plot(x,y,'g', label='Summations')
plt.plot(x,z,'r', label='InitiationRates')
plt.xlabel("Indices")
plt.ylabel("Summations")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
# plt.show()
plt.savefig("Plot of the Summations for " + `csvFile` + ".png")
plt.clf()
#~~Make a plot of Minimum Free Energy~~
y = sortedFeaturesArray[:,-2] #MFE
z = sortedFeaturesArray[:,-1] #Initiation Rates
print actualArray.shape
print y.shape
x = np.arange(len(y))
plt.plot(x,y,'g', label='MFE')
plt.plot(x,z,'r', label='InitiationRates')
plt.xlabel("Indices")
plt.ylabel("MFE")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
# plt.show()
plt.savefig("Plot of the MFE for " + `csvFile` + ".png")
plt.clf()
exit()
'''
# Scale the features
stdsc = preprocessing.StandardScaler()
featureArray_scaled = stdsc.fit_transform(sortedFeaturesArray)

#featureArray_scaled is the feature array scaled and should be used for representation and not for the regression.
print "featureArray_scaled", featureArray_scaled.shape
'''
## Do the colour map ##
myColours = [- HEatmap
            [0,"#ffffff"],   # color of minimun actualArray[:,-1]level = no base
            [0.25,"#0000ff"], [0.5,"#ff0000"], [0.75,"#008000"],  # in-between Blue A, Red T, Green G
            [1,"#ffff00"]   # color of maximum level (from 'zmax'), Yellow C
        ]
forHeatmap = Heatmap(
    z=allNumberedSequences,
    zauto=False, # (!) custom color levels
    zmin=0,      # (!) value of min color level
    zmax=4,      # (!) value of max color levelv
    colorscale=myColours # (!) custom color scales list of lists
    )
data = Data([forHeatmap])
layout = Layout(
    title='Sequences ' + `csvFile`+ ' Whole Length-heatmap'
)
fig = Figure(data=data, layout=layout)
#Make data object
plot_url = py.plot(fig, filename='Sequences ' + `csvFile`+ ' Whole Length-heatmap') #Create it on the browser

## Create the HeatMap ##
#~~~~~~~~~~~~~~
temp = np.array(featureArray_scaled[:,1:5],dtype=float)  # 0=Length, 1,2,3,4 are the frequencies
maxValue=temp.max(axis=0) #axis=0:column: Indicates the no of Points/rows.
minValue=temp.min(axis=0)
print "maxValue for bases freq", maxValue
print "minValue for bases freq", minValue

forHeatmap = Heatmap(
    z=featureArray_scaled[:,1:5],
    zauto=False, # (!) custom color levels
    zmin=minValue,      # (!) value of min color level
    zmax=maxValue,      # (!) value of max color level
    # colorscale=myColours # (!) custom color scales list of lists
    )
data = Data([forHeatmap])
layout = Layout(
    title='Bases Frequency' + `csvFile`
)
fig = Figure(data=data, layout=layout)
plot_url = py.plot(data, filename='Bases Frequency ' + `csvFile`+ ' heatmap')

#~~~~~~~~~~~~~~
temp = np.array(featureArray_scaled[:,5:21],dtype=float) # 5 ....20 are 2mers
maxValue=temp.max(axis=0) #axis=0:column: Indicates the no of Points/rows.
minValue=temp.min(axis=0)
print "maxValue for 2 mers freq", maxValue
print "minValue for 2 mers freq", minValue
forHeatmap = Heatmap(
    z=featureArray_scaled[:,5:21],
    zauto=False, # (!) custom color levels
    zmin=minValue,      # (!) value of min color level
    zmax=maxValue,      # (!) value of max color level
    # colorscale=myColours # (!) custom color scales list of lists
    )
data = Data([forHeatmap])
layout = Layout(
    title='2Mers Frequency ' + `csvFile`
)
fig = Figure(data=data, layout=layout)
plot_url = py.plot(data, filename='2Mers Frequency ' + `csvFile`+ ' heatmap')
#~~~~~~~~~~~~~~
temp = np.array(featureArray_scaled[:,21:85],dtype=float)  # 21....84 are the 3mers
maxValue=temp.max(axis=0) #axis=0:column: Indicates the no of Points/rows.
minValue=temp.min(axis=0)
print "maxValue for 3 mers freq", maxValue
print "minValue for 3 mers freq", minValue
forHeatmap = Heatmap(
    z=featureArray_scaled[:,21:85],
    zauto=False, # (!) custom color levels
    zmin=minValue,      # (!) value of min color level
    zmax=maxValue,      # (!) value of max color level
    # colorscale=myColours # (!) custom color scales list of lists
    )
data = Data([forHeatmap])
layout = Layout(
    title='3Mers Frequency ' + `csvFile`
)
fig = Figure(data=data, layout=layout)
plot_url = py.plot(data, filename='3Mers Frequency ' + `csvFile`+ ' heatmap')
# #~~~~~~~~~~~~~~
exit()
'''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
#THe Init Rates are: featureArray[:,-1]
#Divide the initiation rates into 3 subgroups , 1-2-3 : low rates till -11.5, from -11.5 till -9, from -9 till max.
#Why to do that? To simplify the problem. Do classification instead of Regression.
noOfSeqs = len(featureArray)
for seq in range(0,noOfSeqs):
    if featureArray[seq,-1] < -11.5:
        featureArray[seq,-1] = 1
    else:
        if featureArray[seq,-1] > -9:
            featureArray[seq,-1] = 3
        else:
            featureArray[seq,-1] = 2
    print featureArray[seq][-1]
#This Solution brings too good results, 'cause the 2nd group is way too big.
'''

#Divide the Initiation Rates into classes (many classes)
#Apply the Cross Validation on the different classes
#Use regression to predict the original Initiation Rates.
#Something like that:
#[min,-12][-12,-11.5][-11.5,-11][-11,-10.5][-10.5,-10][-10,-9][-9,-8][-8,6][-6,-4][-4,max]

noOfSeqs = len(featureArray)

# Add an extra column of 0s 
b = np.zeros(noOfSeqs)
instanceArray = np.vstack([featureArray.T,b]).T
featureArray = instanceArray

#Now divide the Initiation Rates into 10 groups
#Now pass from all sequences and replace the last value(0) with a number from 1 to 10
for seq in range(0,noOfSeqs):
    if featureArray[seq,-2] <= np.percentile(featureArray[:,-2],10):
        featureArray[seq,-1] = 1
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],10) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],20):
        featureArray[seq,-1] = 2
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],20) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],30):
        featureArray[seq,-1] = 3
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],30) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],40):
        featureArray[seq,-1] = 4
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],40) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],50):
        featureArray[seq,-1] = 5
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],50) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],60):
        featureArray[seq,-1] = 6
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],60) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],70):
        featureArray[seq,-1] = 7
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],70) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],80):
        featureArray[seq,-1] = 8
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],80) and featureArray[seq,-2] <= np.percentile(featureArray[:,-2],90):
        featureArray[seq,-1] = 9
    elif featureArray[seq,-2] > np.percentile(featureArray[:,-2],90):
        featureArray[seq,-1] = 10


#~~~~~~~~~~~~~~~~~REGRESSION~~~~~~~~~~~~~~~~~
start_time = time.time()
'''
##### Define X and y
X = featureArray[:,0:-1] # features
print "X", X.shape
# print X
y = featureArray[:,-1] # Initiation Rates
print "y", y.shape
# print y
'''
#Do the regression per location
X = []
y = []



for c in range(1,18): #1 to 17 chromosomes.
    #Take all the rows per chromosome
    allcChromosomes = featureArray[featureArray[:,0] == float(c)] # np.where(featureArray[:,1]==0)
    # allcChromosomes = array(chromosome, 87 features, initRates, Class) allcChromosomes(?,90)
    # print "allcChromosomes", allcChromosomes.shape

    X = allcChromosomes[:,1:-2] #Features
    y = allcChromosomes[:,-2] #initiation rates
    # print "X", X.shape
    # print "y", y.shape

    ##### Calculate the Correlation among the features and the target

    pearsonsCorrelations = []
    spearmanCorrelations = []

    #Check if there is at least one sequence.
    if len(y) != 0:
        for i in range(0, len(X[0])): #per column.
            pC = scipy.stats.pearsonr(X[:,i], y)
            pearsonsCorrelations.append(pC)
            # print pC
            sC = scipy.stats.spearmanr(X[:,i], y)
            spearmanCorrelations.append(sC)
            # print sC
        # pCorrelation = scipy.stats.pearsonr(X[:,0], y)
        # print pCorrelation
        # print "pearsonsCorrelations", pearsonsCorrelations
        # print "spearmanCorrelations", spearmanCorrelations

        pearsonsCorrelations = np.array(pearsonsCorrelations)
        spearmanCorrelations = np.array(spearmanCorrelations)
        # print "pearsonsCorrelations", pearsonsCorrelations.shape
        # print "spearmanCorrelations", spearmanCorrelations.shape

        h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
        z = spearmanCorrelations[:,0] # Spearman
        x = np.arange(len(h)) # Features

        plt.figure()
        plt.plot(x,h,'g', label='pearsonsCorrelations')
        plt.plot(x,z,'r', label='spearmanCorrelations')
        plt.xlabel("Features")
        plt.ylabel("Correlations")
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)

        plt.grid(True)
        plt.suptitle('Correlation Coefficient, chromosome:'+ `c`)
        art=[]
        lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
        art.append(lgd)
        # plt.show()
        plt.savefig("ThresholdedIntRates/CorrelationCoefficientFor" + `c`+"Chrom"+ `csvFile` + ".png", additional_artists=art, bbox_inches='tight')
        plt.clf()
        print "------------------------------------------------------------finished Correlations Calculation"
        #~~~~~~~~~~Cross Valication ~~~~~~~~~~~~~~~~~~~~~~~~~~
        skf = cross_validation.KFold(len(y),n_folds=3)
        #skf = cross_validation.KFold(allcChromosomes[:,-1],n_folds=5)
        #skf = cross_validation.StratifiedKFold(allcChromosomes[:,-1], n_folds=3) #Stratified Cross Validation


        # Do SVR
        allScores = []
        print "X and Y", X.shape, y.shape
        fold = 1
        for train_index, test_index in skf:
            #print("TRAIN:", train_index, "TEST:", test_index)
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            # svc_linear = svm.SVC(kernel='linear', degree=4, C=1, cache_size=300)
            svr_linear = svm.SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma=0.0,
                kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)        #SupportVectorMachine - Support Vector Regression
            print X_train.shape
            print y_train.shape

            start_training_time = time.time()
            print "start training"
            svr_linear.fit(X_train, y_train)
            print "Done training in ", time.time() - start_training_time

            score = svr_linear.score(X_test, y_test)
            allScores.append(score)
            y_linear_predictTest = svr_linear.predict(X_test)
            y_linear_predictTrain = svr_linear.predict(X_train)

            ys = np.vstack([y_train , y_linear_predictTrain]).T
            ys = ys[ys[:, 0].argsort()]
            #print ys
            y_train = ys[:,0]
            y_linear_predictTrain = ys[:,1]

            ys = np.vstack([y_test , y_linear_predictTest]).T
            ys = ys[ys[:, 0].argsort()]
            #print ys
            y_test = ys[:,0]
            y_linear_predictTest = ys[:,1]

            #Images
            plt.close('all')
            plt.figure()
            plt.subplots_adjust(hspace=.5, wspace=.5)
            plt.subplot(1,2,1)
            plt.title('SVR train')
            allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=3, label='train')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_linear_predictTrain, c='r', edgecolor='none', s=3, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('Init Rates 1,2,3')#('initiation rates (log)')
            art=[]
            lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
            art.append(lgd)

            plt.subplot(1,2,2)
            allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=3, label='test')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_linear_predictTest, c='r', edgecolor='none', s=3, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('Init Rates 1,2,3')#('initiation rates (log)')
            plt.title('SVR test')
            lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
            art.append(lgd)
            plt.suptitle('SVR_AllFeats'+ `c` + 'Chrom'+`fold`+'Fold')
            plt.savefig('ThresholdedIntRates/rbf1_AllFeatsFor' + `c` + 'Chrom' + `csvFile` +`fold`+'.png', additional_artists=art, bbox_inches='tight')
            # plt.show()
            plt.clf()
            fold += 1



        # Print SVR scores
        allScores = np.array(allScores)
        print "SVR Scores:" + `c`, allScores
        print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
        print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_linear_predictTrain, y_train)))
        print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_linear_predictTest, y_test)))


        #Do Random Forest
        allScores = []
        fold = 1
        for train_index, test_index in skf:
            #print("TRAIN:", train_index, "TEST:", test_index)
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            rndmForest = RandomForestRegressor(n_estimators=100)
            # rndmForest = RandomForestClassifier(n_estimators=100)
            print "Start training Random Forest"
            start_time_training = time.time()
            rndmForest.fit(X_train, y_train)
            print "Done training in ",time.time() - start_time_training, "seconds"
            score = rndmForest.score(X_test, y_test)
            allScores.append(score)
            y_predictTest = rndmForest.predict(X_test)
            y_predictTrain = rndmForest.predict(X_train)

            ys = np.vstack([y_train , y_predictTrain]).T
            ys = ys[ys[:, 0].argsort()]
            #print ys
            y_train = ys[:,0]
            y_predictTrain = ys[:,1]

            ys = np.vstack([y_test , y_predictTest]).T
            ys = ys[ys[:, 0].argsort()]
            #print ys
            y_test = ys[:,0]
            y_predictTest = ys[:,1]

            #
            plt.close('all')
            plt.figure()
            plt.subplots_adjust(hspace=.5, wspace=.5)
            plt.subplot(1,2,1)
            plt.title('rndmForest train scaled')
            allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=3, label='train')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_predictTrain, c='r', edgecolor='none', s=3, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')
            art = []
            lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
            art.append(lgd)

            plt.subplot(1,2,2)
            allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=3, label='test')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_predictTest, c='r', edgecolor='none', s=3, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')
            plt.title('rndmForest test')
            lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
            art.append(lgd)
            plt.suptitle('rndmForest_AllFeats'+ `c` + 'Chrom'+`fold`+'Fold')
            plt.savefig('ThresholdedIntRates/rndmForest_AllFeatsFor' + `c` + 'Chrom' + `csvFile` +`fold`+'.png', additional_artists=art, bbox_inches='tight')
            # plt.show()
            plt.clf()
            fold += 1

        # print "Prediction of X is ", Y_predicted
        allScores = np.array(allScores)
        print "Random Forest Scores:" + `c`, allScores
        print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
        print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_predictTrain, y_train)))
        print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_predictTest, y_test)))
        print "Feature Importance for Chromosome: "+ `c`, rndmForest.feature_importances_

        importances = rndmForest.feature_importances_
        std = np.std([rndmForest.feature_importances_ for rndmForest in rndmForest.estimators_],
             axis=0)
        # print "std",std
        indices = np.argsort(importances)[::-1]
        # print "indices",indices

        # print "Feature Ranking For Chromosome: " + `c`
        # for f in range(1,88): # I have 87 features
        #     print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

        # Plot the feature importances of the forest
        plt.figure(figsize=(26, 12))
        plt.title("Feature importances")
        plt.bar(range(1,88), importances,
               color="r", yerr=std, align="center")
        plt.xticks(np.arange(1,88))
        plt.xlim([-1, 88])
        plt.grid(True)
        plt.savefig("ThresholdedIntRates/FeatImportChrom"+`c`+".png")
        #plt.show()




print "Regression takes"
print time.time() - start_time



exit()

