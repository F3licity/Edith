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
##########################################################################################################################
#######Chromosomes

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
# print "chromosomes", chromosomes[0]
# exit()
#######################################################################
#When searching pairs of base pairs or 3 base pairs
#It can easily be used for more, by adding a for
def find_kmers(sequence):

    #Create the k mers arrays
    dimers = []
    trimers = []
    bases = ["A","T","G","C"]
    for i in bases:
        for j in bases:
            dimers.append(i+j)
            for y in bases:
                trimers.append(i+j+y)

    #Count the occurences of the kmers in the sequences
    dimercounts = []
    for i in dimers:
        dimercounts.append(sequence.count(i))

    trimercounts =[]
    for j in trimers:
        trimercounts.append(sequence.count(j))

    return dimers,trimers,dimercounts,trimercounts

#because the chromosomes are numbered with I, II etc for my own ease I use 1,2 etc
def romanToNumeric(source):
    #Roman to Numerical
    if source == 'chrI':
        source = 1
    elif source == 'chrII':
        source = 2
    elif source == 'chrIII':
        source = 3
    elif source == 'chrIV':
        source = 4
    elif source == 'chrV':
        source = 5
    elif source == 'chrVI':
        source = 6
    elif source == 'chrVII':
        source = 7
    elif source == 'chrVIII':
        source = 8
    elif source == 'chrIX':
        source = 9
    elif source == 'chrX':
        source = 10
    elif source == 'chrXI':
        source = 11
    elif source == 'chrXII':
        source = 12
    elif source == 'chrXIII':
        source = 13
    elif source == 'chrXIV':
        source = 14
    elif source == 'chrXV':
        source = 15
    elif source == 'chrXVI':
        source = 16
    elif source == 'chrmt':
        source = 17

    return source

# The following function searches a given term inside an array in a specific column and 
#returns True if found, or False otherwise.
def searchIfInColumn(term, array, column):

    for i in range(len(array)):
        if term == array[i][column]:
            return True
    return False

# The following function searches a given term inside an array in a specific column and 
#returns the Index if found, or -1 otherwise.
def searchIndex(term, array, column, entry):

    for i in range(len(array)):
        if entry == '':
            if term == array[i][column]:
                return i
        else:
            if term in array[i][column][entry]:
                return i
    return -1
###################################################################################################

#Open also the yeast names file withcsvFile = 'newDataset.csv' the translation initiation rates
csvFile = 'newDataset.csv'
# csvFile = 'ciandrini.csv'
'''
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

genes           = np.array(genes)
InitiationRates = np.array(InitiationRates)
print "genes: ",           genes.shape
print "InitiationRates: ", InitiationRates.shape
# print genes[0]
# print InitiationRates
print "--------------------------------------------------------------------------------------------------end of csv file."
'''
########################################################################
'''
# Import GFF3 Files
nagalakshmiData = np.load("nagalakshmiData.npy")
nagalakshmiData = np.array(nagalakshmiData)
print "nagalakshmiData", nagalakshmiData.shape

yassourData = np.load("yassourData.npy")
yassourData = np.array(yassourData)
print "yassourData", yassourData.shape

union = np.vstack((nagalakshmiData,yassourData))
# print "union[0]", union[0].T.shape
print "InitiationRates,shape", InitiationRates[0].shape
# print union[0][0]
# exit()
'''
#######################################################################
# Save in preArray all the info when matching Genes and Gff3 Files.
'''
preArray = []
for g in genes:
    for u in range(len(union)):
        if g in union[u][8]['Name'][0:7]:
            geneIndex = np.where(genes == g) #type: tuple
            preArray.append(np.hstack((union[u],InitiationRates[geneIndex[0][0]])))
np.save("preArray", preArray)
'''

preArray = np.load("preArray.npy")
print "preArray",preArray.shape
# print "preArray", preArray[0:5]
# print preArray[16][8]['Name']
# print preArray[20][8]['Name']


###########################################################################
'''
#Create an array for the features
#Rows as many as the genes that will be used
#Columns 5: 5th_start, 5th_end, 3rd_start, 3rd_end

finalGenes = []
duplicates = []
maxDuplicateOccurence = 0;
#finalGenes : Name, 5start, 5end, 3start, 3end, chromosome, TransInitRates
for i in range(len(preArray)):

    #If it is not already in the array, add it
    if (searchIfInColumn(preArray[i][8]['Name'][:-5], finalGenes, 0) is False):

        #check whether it is 5 UTR or 3 UTR
        if preArray[i][8]['Name'].find('5') !=-1: #if its 5utr
            gene = [preArray[i][8]['Name'][:-5], preArray[i][3], preArray[i][4], 0, 0,preArray[i][0],0]
            finalGenes.append(gene)

        else: #if its 3utr
            gene = [preArray[i][8]['Name'][:-5], 0, 0, preArray[i][3], preArray[i][4],preArray[i][0],0]
            finalGenes.append(gene)

    #If it is in the array, update the missing values
    else:
        # find gene's position in the array.
        geneIndex = searchIndex(preArray[i][8]['Name'][:-5],finalGenes,0,'')


        # if what we have already is 3utr save the 5utr info
        if (finalGenes[geneIndex][1] == 0) and (preArray[i][8]['Name'].find('5UTR') !=-1):
            # if the 5utr is smaller than the 3utr then, complte the info
            if preArray[i][4] < finalGenes[geneIndex][3] :
                finalGenes[geneIndex][1] = preArray[i][3]
                finalGenes[geneIndex][2] = preArray[i][4]

        # if what we have already is 5utr save the 3utr info
        elif (finalGenes[geneIndex][3] == 0) and (preArray[i][8]['Name'].find('3UTR') !=-1):
            # if the 3utr start is bigger than the 5utr end, complete the info
            if preArray[i][3] > finalGenes[geneIndex][2] :
                finalGenes[geneIndex][3] = preArray[i][3]
                finalGenes[geneIndex][4] = preArray[i][4]

        else: # see duplicates
            duplicates.append(preArray[i][8]['Name'][:-5])
            duplicateCounter = duplicates.count(preArray[i][8]['Name'][:-5])
            if (duplicateCounter > maxDuplicateOccurence):
                maxDuplicateOccurence = duplicateCounter

            # if you find 5_utr_end smaller, update
            if finalGenes[geneIndex][2] > preArray[i][4]:
                finalGenes[geneIndex][1] = preArray[i][3]
                finalGenes[geneIndex][2] = preArray[i][4]
            #if you find 3Utr start bigger, update.
            if finalGenes[geneIndex][3] < preArray[i][3]:
                finalGenes[geneIndex][3] = preArray[i][3]
                finalGenes[geneIndex][4] = preArray[i][4]


finalGenes = np.array(finalGenes)
print "finalGenes", finalGenes.shape
# print "finalGenes", finalGenes[0:10]
# print "finalGenes",finalGenes

duplicates = np.array(duplicates)
#print "dupliacets:",duplicates #np.unique(duplicates)
print "duplicates", np.unique(duplicates).shape

print maxDuplicateOccurence

#Add the Initiation Rates as the last column
for i in range(len(finalGenes)):

    pos = finalGenes[i][0].find('-')
    if pos != -1:
        name = finalGenes[i][0][:pos]
    else: # if -1 means there is not
        name = finalGenes[i][0]

    geneIndex =  np.where(genes==name)[0]

    if geneIndex != -1:
        finalGenes[i][-1] = InitiationRates[geneIndex][0]

# print "finalGenes", finalGenes[0:15]
print "finalGenes",finalGenes.shape


#Now split the array, into Name and Features
arrayNames = finalGenes[:,0]  # first column
arrayFeatures = finalGenes[:,1:] # all the rest columns

# print "arrayNames.shape",arrayNames[0]
# print "arrayFeatures.shape",arrayFeatures[0]

# print "arrayFeatures", arrayFeatures[0]
# print np.nonzero(arrayFeatures[0])

# Get the Length of the Coding Sequences
arrayCodSequences = []
codingsLengths = []
count = 0
for i in range(len(arrayNames)):

    #Roman to Numerical to find which Chromosome
    chromosomeNo = romanToNumeric(arrayFeatures[i][4])
    # print "chromosomeNo",chromosomeNo

    #Coding Sequence : 5utr end  till 3 utr start
    codingStart = int(arrayFeatures[i][1])
    codingEnd   = int(arrayFeatures[i][2])

    # codingStart = int(arrayFeatures[i][0])
    # codingEnd = int(arrayFeatures[i][3])

    # print "codingStart", codingStart
    # print "codingEnd", codingEnd
    if codingEnd < codingStart:
        count +=1
        codingSequence = ''
    else: 
        codingSequence = chromosomes[chromosomeNo][codingStart:codingEnd]
        if (len(codingSequence )==0):
            print "why? ",codingStart,chromosomeNo,codingEnd, len(chromosomes[chromosomeNo])
    arrayCodSequences.append(codingSequence)
    # print "codingSequence", codingSequence
    # print "len(codingSequence) "+`i`+" : ", len(codingSequence)
    codingsLengths.append(len(codingSequence))


print "Sequences that have codingEnd < codingStart are: ", count


codingsLengths    = np.array(codingsLengths)
arrayCodSequences = np.array(arrayCodSequences)
print "arrayNames : ", arrayNames.shape
print "codingsLengths : ", codingsLengths.shape
print "arrayCodSequences : ",arrayCodSequences.shape

# # print "arrayNames", arrayNames[0:8]
# print "arrayFeatures", arrayFeatures[0:15]
# print "arrayCodSequences[0:15]" , arrayCodSequences[0:15]
# # for n in range(len(arrayCodSequences)):
# #     print arrayCodSequences[n][0:10]
# print "codingsLengths", codingsLengths[0:15]
'''

'''
mask = np.ones(len(arrayNames), dtype=bool)
for i in range(len(arrayNames)):

    #If there is no 5utr info discard Gene
    if arrayFeatures[i][0] == '0':
        mask[i] = False

    #If the length of the gene is 0, discard this Gene
    if codingsLengths[i] == 0:
        mask[i] = False


arrayFeatures = arrayFeatures[mask]
arrayNames = arrayNames[mask]
arrayCodSequences = arrayCodSequences[mask]
codingsLengths = codingsLengths[mask]

np.save("arrayNames", arrayNames)
np.save("arrayFeatures", arrayFeatures)
np.save("codingsLengths", codingsLengths)
np.save("arrayCodSequences", arrayCodSequences)
'''
arrayNames        = np.load("arrayNames.npy")
arrayFeatures     = np.load("arrayFeatures.npy")
codingsLengths     = np.load("codingsLengths.npy")
arrayCodSequences = np.load("arrayCodSequences.npy")

print "arrayNames : ", arrayNames.shape
print "arrayFeatures : ", arrayFeatures.shape
print "codingsLengths : ", codingsLengths.shape
print "arrayCodSequences : ",arrayCodSequences.shape

X = codingsLengths.astype(float)
y = arrayFeatures[:,-1].astype(float)

pCorrelation = scipy.stats.pearsonr(X,y)
sCorrelation = scipy.stats.spearmanr(X, y)
print "Correlation among Coding's Length and Init Rates"
print "Pearson : ", pCorrelation
print "Spearman : ", sCorrelation
exit()
skf = cross_validation.KFold(len(y),n_folds=5)

# # Do SVR
# allScores = []
# print "X and Y", X.shape, y.shape
# fold = 1
# for train_index, test_index in skf:
#     #print("TRAIN:", train_index, "TEST:", test_index)
#     X_train, X_test = X[train_index], X[test_index]
#     y_train, y_test = y[train_index], y[test_index]
#     print "X_train.shape",X_train.shape
#     print "y_train.shape", y_train.shape
#     svr_linear = svm.SVR(kernel='linear', degree=4, C=1, cache_size=300)
#     # svr_linear = svm.SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma=0.0,
#     #     kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)        #SupportVectorMachine - Support Vector Regression
#     print "start training"
#     X_train = X_train.reshape( (len(X_train), 1) )
#     svr_linear.fit(X_train, y_train)
#     print "done training"
#     score = svr_linear.score(X_test, y_test)
#     allScores.append(score)
#     y_linear_predictTest = svr_linear.predict(X_test)
#     y_linear_predictTrain = svr_linear.predict(X_train)

#     ys = np.vstack([y_train , y_linear_predictTrain]).T
#     ys = ys[ys[:, 0].argsort()]
#     #print ys
#     y_train = ys[:,0]
#     y_linear_predictTrain = ys[:,1]

#     ys = np.vstack([y_test , y_linear_predictTest]).T
#     ys = ys[ys[:, 0].argsort()]
#     #print ys
#     y_test = ys[:,0]
#     y_linear_predictTest = ys[:,1]

#     #Images
#     plt.close('all')
#     plt.figure()
#     plt.subplots_adjust(hspace=.5, wspace=.5)
#     plt.subplot(1,2,1)
#     plt.title('SVR train')
#     allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
#     plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=3, label='train')
#     plt.hold('on')
#     plt.scatter(allSampleIndexes, y_linear_predictTrain, c='r', edgecolor='none', s=3, label='prediction')
#     plt.xlabel('#of rows')
#     plt.ylabel('initiation rates (log)')
#     art=[]
#     lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
#     art.append(lgd)

#     plt.subplot(1,2,2)
#     allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
#     plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=3, label='test')
#     plt.hold('on')
#     plt.scatter(allSampleIndexes, y_linear_predictTest, c='r', edgecolor='none', s=3, label='prediction')
#     plt.xlabel('#of rows')
#     plt.ylabel('initiation rates (log)')
#     plt.title('SVR test')
#     lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
#     art.append(lgd)
#     plt.suptitle('SVR_length'+ `fold`+'Fold')
#     plt.savefig('CodingRegression/svr_length'+ `csvFile` +'Fold '+`fold`+'.png', additional_artists=art, bbox_inches='tight')
#     # plt.show()
#     plt.clf()
#     fold += 1



# # Print SVR scores
# allScores = np.array(allScores)
# print "SVR Scores:", allScores
# print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
# print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_linear_predictTrain, y_train)))
# print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_linear_predictTest, y_test)))

#Do Random Forest
allScores = []
fold = 1
for train_index, test_index in skf:
    #print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    X_train = X_train.reshape( (len(X_train), 1) )
    X_test  = X_test.reshape((len(X_test), 1))
    rndmForest = RandomForestRegressor(n_estimators=100)
    print "Start training Random Forest"
    rndmForest.fit(X_train, y_train)
    print "Done training"
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

    #Images
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
    plt.suptitle('rndmForest_length'+`fold`+'Fold')
    plt.savefig('CodingRegression/rndmForest_length'+ `csvFile` +'Fold '+`fold`+'.png', additional_artists=art, bbox_inches='tight')
    # plt.show()
    plt.clf()
    fold += 1

# print "Prediction of X is ", Y_predicted
allScores = np.array(allScores)
print "Random Forest Scores:", allScores
print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_predictTrain, y_train)))
print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_predictTest, y_test)))





