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

#######################################################################
#######Chromosomes

chromosomes = []

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

#Now I have an array of Chromosomes that has all the 17 chromosomes 0....16
#print chromosomes[16]
print "--------------------------------------------------------------------------------------------------end of chromosomes"
#######################################################################

#Open also the yeast names file withcsvFile = 'newDataset.csv' the translation initiation rates
csvFile = 'newDataset.csv'
# csvFile = 'ciandrini.csv'
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
# print genes
# print InitiationRates
print "--------------------------------------------------------------------------------------------------end of csv file."
########################################################################
#Parse gff3 file.



#Initialized GeneInfo named tuple. Note: namedtuple is immutable(ametavlhtos)
gffInfoFields = ["source", "filename", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

#------------------------------------------------------------------------------------------
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dictionary"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.unquote(key)] = urllib.unquote(value)
    return ret

def parseGFF3(filename):

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

if __name__ == "__main__":
    nagalakshmiData = []
    yassourData = []

    for record in parseGFF3("Nagalakshmi_2008_UTRs.gff3"):
        nagalakshmiData.append(record)
    for record in parseGFF3("Yassour_2009_UTRs.gff3"):
        yassourData.append(record)

    nagalakshmiData = np.array(nagalakshmiData)
    yassourData = np.array(yassourData)

    print "nagalakshmiData", nagalakshmiData.shape
    print nagalakshmiData[:,8]
    print "yassourData", yassourData.shape
    print yassourData[0]

    np.save("nagalakshmiData", nagalakshmiData)
    np.save("yassourData", yassourData)
    exit()


 













    recordCount = 0
    five  = []
    three = []
    

    for record in parseGFF3("Nagalakshmi_2008_UTRs.gff3"):
    	wholeName = record.attributes['Name']

    	if '5UTR' in wholeName:
    		five.append(wholeName[:-5])
    	elif '3UTR' in wholeName:
    		three.append(wholeName[:-5])
    	else:
    		print "record", record

    five  = np.array(five)
    # print "five", five
    three = np.array(three)
    # print "three", three
    print "five", five.shape
    print "three", three.shape

    both = list(set(five) & set(three)) #Intersection
    # print "both", both
    both = np.array(both)
    print "both", both.shape
    

    notInFive = list(set(three).difference(set(five)))
    print "notInFive", notInFive

    notInThree = list(set(five).difference(set(three)))
    print "notInThree", notInThree
    exit()

    differences = list(set(three).symmetric_difference(set(five)))
    print "differences", differences
    differences = np.array(differences)
    print "differences",differences.shape

    exceptions = []
    for i in three:
        pos = i.find('-')

        if pos != -1:
            i = i[:pos]
        else:
            i = i

        if i in genes:
            exceptions.append(i)

    print "exceptions", exceptions
    exceptions = np.array(exceptions)
    print "exceptions shape : ", exceptions.shape

    exit()

    extraThrees = []
    for t in three:
        matching = [s for s in five if t[:-5] in s]
        print "matching",matching
        
    # 	if t[:-5] not in five:
    # 		print t[:-5]
    # 		extraThrees.append(t)
    # extraThrees = np.array(extraThrees)
    # print "extraThrees", extraThrees.shape
    exit()
    extraFives = []
    for f in five:
    	if f[:-5] not in three:
    		print f[:-5]
    		extraFives.append(f)
    extraFives = np.array(extraFives)
    print "extraFives", extraFives.shape

    # extragenes = []
    # for f in extraFives:
    # 	if extraFives[f][0:7] not in genes:
    # 		print f, extraFives[f]
    # 		extragenes.append(extraFives[f])
    # print extragenes
    # for t in extraThrees:
    # 	if extraThrees[t][-] not in genes:
    # 		print t, extraThrees[f]
    # 		extragenes.append(extraThrees[f])
    # extragenes = np.array(extragenes)
    # print "extragenes", extragenes.shape