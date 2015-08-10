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
from sklearn.pipeline import Pipeline, make_pipeline
#------------------------------------------------------------------------------------------
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dictionary"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.unquote(key)] = urllib.unquote(value)
    return ret

def parseGFF3(filename, gffInfoFields, GFFRecord):

    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            #print line

            parts = line.strip().split("\t")

            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                #"seqid": None if parts[0] == "." else urllib.unquote(parts[0]), #cause i dont have seqid-
                "source": None if parts[0] == "." else urllib.unquote(parts[0]),
                "filename": None if parts[1] == "." else urllib.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo) #**kwargs allow you to pass a variable number of arguments to a function. 
            #explained in :https://freepythontips.wordpress.com/2013/08/04/args-and-kwargs-in-python-explained/

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

def findBasesFrequency(sequence):
    from collections import Counter

    baseCounter = Counter(sequence)
    sequenceLength = len(sequence)

    Afrequency = baseCounter['A']/float(sequenceLength)
    Tfrequency = baseCounter['T']/float(sequenceLength)
    Gfrequency = baseCounter['G']/float(sequenceLength)
    Cfrequency = baseCounter['C']/float(sequenceLength)

    return Afrequency, Tfrequency, Gfrequency, Cfrequency

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

#I need the first 100bp upstream and not the whole 5UTR
#  ------------------------|--coding sequence---
#                ^^^^^^^^^^
def shortenSequence(start,end):

    newStart = end -99
    if newStart > start:
        start = newStart

    return start

#Find Duplicates in a list
def findDuplcicates(a):
    duplicates = []
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            if a[i] == a[j]:
                duplicates.append(a[i])
    return duplicates


# The following function searches a given term inside an array in a specific column and 
#returns True if found, or False otherwise.
def searchIfInColumn(term, array, column):

    for i in range(len(array)):
        if term == array[i][column]:
            return True
    return False


def plotWith4Subplots(howmany, x1, name1, x2, name2, x3, name3, x4, name4, barWidth, sharex, title, where):
    plt.close('all')
    f, (ax1, ax2, ax3, ax4) = plt.subplots(howmany)
    f.tight_layout()
    ax1.bar(barWidth, x1,0.8,color='b',linewidth=0)
    ax1.set_title(name1)
    ax2.bar(barWidth, x2,0.8,color='r',linewidth=0)
    ax2.set_title(name2)
    ax3.bar(barWidth, x3,0.8,color='g',linewidth=0)
    ax3.set_title(name3)
    ax4.bar(barWidth, x4,0.8,color='y',linewidth=0)
    ax4.set_title(name4)

    #f.xlabel('Direction: 3 to 5')
    f.suptitle(title, y=1.08)
    #plt.show()
    f.savefig(where)
    plt.clf()


#For 2 different measurements
def plotTwoScales(x,x_label, y1,y1_label, y2, y2_label, labelone, labeltwo, filename):

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt

    host = host_subplot(111, axes_class=AA.Axes)
    # plt.subplots_adjust(right=0.75)  # Is for an axis that is a bit further than the plot

    par1 = host.twinx()
    #par2 = host.twinx()

    # offset = 60
    # new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    #par2.axis["right"] = new_fixed_axis(loc="right",
                                        # axes=par2,
                                        # offset=(offset, 0))

    # par2.axis["right"].toggle(all=True)

    # host.set_xlim(0, 2)  #Set limits
    # host.set_ylim(0, 2)

    host.set_xlabel(x_label)
    host.set_ylabel(y1_label)
    par1.set_ylabel(y2_label)
    #par2.set_ylabel("Velocity")

    p1, = host.plot(x, y1, 'ro', label = y1_label)
    p2, = par1.plot(x, y2, 'go', label = y2_label)
    #p3, = par2.plot([0, 1, 2], [50, 30, 15], label="Velocity")

    # par1.set_ylim(0, 4)
    # par2.set_ylim(1, 65)

    host.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                ncol=2, mode="expand", borderaxespad=0.)

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    #par2.axis["right"].label.set_color(p3.get_color())

    plt.draw()
    # plt.show()

    plt.savefig(filename)


#Slide window
class slidingWindow(object):
    '''Returns iterator that will emit chunks of size 'winSize' each time self.next()
    is called.'''
    def __init__(self,sequence,winSize,step=1):
        '''Returns iterator that will emit chunks of size 'winSize' and 'step' forward in
        the seq each time self.next() is called.'''

        # verification code
        # if not type(sequence) == type(''):
        #     raise Exception("**ERROR** type(sequence) must be str.")
        if not ((type(winSize) == type(0)) and (type(step) == type(0))):
            raise Exception("**ERROR** type(winSize) and type(step) must be int.")
        if step > winSize:
            raise Exception("**ERROR** step must not be larger than winSize.")
        if winSize > len(sequence):
            raise Exception("**ERROR** winSize must not be larger than sequence length.")
        self._seq = sequence
        self._step = step
        self._start = 0
        self._stop = winSize

        # # Pre-compute number of chunks to emit
        # numOfChunks = ((len(sequence)-winSize)/step)+1

        # # Do the work
        # for i in range(0,numOfChunks*step,step):
        #     yield sequence[i:i+winSize]

    def __iter__(self):
        return self

    def next(self):
        """Returns next window chunk or ends iteration if the sequence has ended."""
        try:
            assert self._stop <= len(self._seq), "Not True!"
            chunk = self._seq[self._start:self._stop]
            self._start += self._step
            self._stop  += self._step
            return chunk
        except AssertionError:
            raise StopIteration



def doRandmFrstRegression(X, y, skf, saveImageName, correlationImageName): #OR LAsso
    from sklearn.ensemble import GradientBoostingRegressor

    #Do Random Forest
    print "Regression Start"
    print ""
    allScores = []
    fold = 1
    for train_index, test_index in skf:
        #print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        #REgression and Feature Selection
        # clf = linear_model.Lasso(alpha=0.0001, normalize=True)
        clf = GradientBoostingRegressor(n_estimators=500, max_depth=4, learning_rate=0.1, loss='huber', min_samples_leaf=3,random_state=0, verbose=1)


        # clf = RandomForestRegressor(n_estimators=1000)
        #rndmForest = make_pipeline(RandomForestRegressor(n_estimators=1000), clf)
        print "X_train",X_train.shape
        print "y_train",y_train.shape

        print "Start training Regressor"
        # start_time_training = time.time()
        clf.fit(X_train, y_train)
        # print "Done training in ",time.time() - start_time_training, "seconds"

        print ("X_test:", X_test.shape, "y_test:", y_test.shape)
        score = clf.score(X_test, y_test)
        print "score: ", score
        allScores.append(score)

        y_predictTest = clf.predict(X_test)
        y_predictTrain = clf.predict(X_train)

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
        plt.title('train set')
        allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
        plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=3, label='train')
        plt.hold('on')
        plt.scatter(allSampleIndexes, y_predictTrain, c='r', edgecolor='none', s=3, label='prediction')
        plt.xlabel('# of sequences')
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
        plt.title('test set')
        lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
        art.append(lgd)
        plt.suptitle('Regression'+`fold`+'Fold')
        plt.savefig(saveImageName +'Fold '+`fold`+'.png', additional_artists=art, bbox_inches='tight')
        # plt.show()
        plt.clf()

        # #Calculate the Correlation at each fold
        # pearsonsCorrelations = []
        # spearmanCorrelations = []
        
        # for i in range(len(X_train.T)): #per column.
        #     #Print only those that are taken into account and not all of them
        #     if clf.coef_[i] != 0:
        #         print "important features" , i
        #     pC = scipy.stats.pearsonr(X_train[:,i], y_train)
        #     pearsonsCorrelations.append(pC)
        #     # print pC
        #     sC = scipy.stats.spearmanr(X_train[:,i], y_train)
        #     spearmanCorrelations.append(sC)
        #     # print sC

        #     pearsonsCorrelations = np.array(pearsonsCorrelations)
        #     spearmanCorrelations = np.array(spearmanCorrelations)
        #     # print "pearsonsCorrelations", pearsonsCorrelations.shape
        #     # print "spearmanCorrelations", spearmanCorrelations.shape

        #     h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
        #     z = spearmanCorrelations[:,0] # Spearman
        #     x = np.arange(len(X_train.T)) # Features
        #     p = pearsonsCorrelations[:,1] # Pearson P values
        #     q = spearmanCorrelations[:,1] # Spearman P values

        #     plt.clf()
        #     ind = np.arange(len(X_train.T)) #width of a bar
        #     f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
        #     f.tight_layout()
        #     ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
        #     ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
        #     ax1.set_xticks(x)
        #     ax1.grid()
        #     ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
        #     ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
        #     ax2.set_title('P value')
        #     ax2.set_xticks(x)
        #     ax2.grid()
        #     f.suptitle('Correlation at fold '+ `fold`)
        #     f.savefig(correlationImageName + `fold`+'.svg')
        #     plt.clf()

        # print("clf.coef_", clf.coef_)
        # print("intercept_", clf.intercept_)
        print "End of Fold", fold
        fold += 1


    # print "Prediction of X is ", Y_predicted
    allScores = np.array(allScores)
    print "Gradient Boosting Scores:", allScores
    print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
    print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_train, y_predictTrain)))
    print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_test, y_predictTest)))


def readCSV(csvFile, delimiter):
    a = []
    with open(csvFile,'rb') as csvfile: #read binary file
        data = list(csv.reader(csvfile, delimiter=delimiter))
        skiprow = False
        
        teller = 0
        for row in data:
            teller += 1
            if len(row) > 1:
                #print ', '.join(row)
                for i in range(len(row)):
                    # if "-" in row[0]:  # That is for reading Yeast Names
                    #     row[0] = row[0][:-2]
                    if "#" in row[i]:
                        skiprow = True
                    #print row
                if skiprow == False:
                    a.append(row)
            else:
                print "row; ",teller,row
    return a

'''
def GetStressValues():

for i in range(noRows): # no or Rows len(actualArrayWithNames) #3409
    
    #IF Gene has Stress Value 
    stressIndexes = np.where(stressedGenes[:,0]==actualArrayWithNames[i][0])
    countValues = len(stressIndexes[0])
    
    #Add as Feature the No of Stress Terms that each gene has
    actualArrayWithNames[i][-2] = countValues  #Here is Stress added as feature
print "Stress added as Feature actualArrayWithNames", actualArrayWithNames.shape

temp = np.array(actualArrayWithNames[:,-2],dtype=int)
maxGroups=temp.max(axis=0)
print ""
print "--maxGroups of Stress that a Gene belongs", maxGroups
print np.where(actualArrayWithNames[:,-2]==maxGroups)
minGroups=temp.min(axis=0)
print "--minGroups of Stress that a Gene belongs", minGroups
print ""
'''
# def calculateCorrelation(X, y):

#     pearsonsCorrelations = []
#     spearmanCorrelations = []

#     for i in range(len(X.T)): #per column.

#         pC = scipy.stats.pearsonr(X[:,i], y)
#         pearsonsCorrelations.append(pC)
#         # print pC
#         sC = scipy.stats.spearmanr(X[:,i], y)
#         spearmanCorrelations.append(sC)
#         # print sC

#         pearsonsCorrelations = np.array(pearsonsCorrelations)
#         spearmanCorrelations = np.array(spearmanCorrelations)
#         # print "pearsonsCorrelations", pearsonsCorrelations.shape
#         # print "spearmanCorrelations", spearmanCorrelations.shape

#         # #Select the interesting ones
#         # for j in range(len(X.T)):
            
#         #     h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
#         #     z = spearmanCorrelations[:,0] # Spearman
#         #     x = np.arange(len(X.T)) # Features
#         #     p = pearsonsCorrelations[:,1] # Pearson P values
#         #     q = spearmanCorrelations[:,1] # Spearman P values

#         #     plt.clf()
#         #     ind = np.arange(len(X.T)) #width of a bar
#         #     f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
#         #     f.tight_layout()
#         #     ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
#         #     ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
#         #     ax1.set_xticks(x)
#         #     ax1.grid()
#         #     ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
#         #     ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
#         #     ax2.set_title('P value')
#         #     ax2.set_xticks(x)
#         #     ax2.grid()
#         #     f.suptitle('Correlation and P Value')
#         #     f.savefig(correlationImageName +'.svg')
#         #     plt.clf()