#In this file I calculate the Conditional Entropy per base and per 3bases.
#As well as the SUmmation of the values that I get from the first.

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
from myfunctions import parseGFFAttributes, parseGFF3, find_kmers, romanToNumeric, shortenSequence, findDuplcicates, searchIfInColumn, plotWith4Subplots, plotTwoScales
#############################################################################################################################################
actualArrayWithNames = np.load('Ciandrini/actualArrayWithNames.npy')
print "actualArrayWithNames Loaded", actualArrayWithNames.shape

sortedArray = []
initRates= actualArrayWithNames[:,-1].astype(float)
sortedArray = actualArrayWithNames[initRates.argsort()] # Sort the array by InitiationRate
sortedArray = sortedArray[:,1:]


#sortedArray = np.load('Gritsenko/sortedArray.npy')
print "sortedArray  loaded, is the actualArrayWithNames Array (without the names) sorted by Init Rates: ", sortedArray.shape

temp = [] # Gather the Lengths (100upstream + 40CDS)
for i in range(len(sortedArray)):
	temp.append(len(sortedArray[i][0]))
temp = np.array(temp)
maxLength=temp.max(axis=0)
print "maxLength", maxLength # Find max length
minLength=temp.min(axis=0)
print "minLength", minLength
print ""

'''
arrayEsetA = []
arrayEsetT = []
arrayEsetG = []
arrayEsetC = []

#Pass from every position
for k in range(maxLength):
	#print "position", k
	tempSequences = []

	#Select the sequences that have this length
	for s in range(len(sortedArray)):

		#Find Sequence's Length
		lenseq = len(sortedArray[s][0]) #int(sortedArray[s][2])
		#print "lenseq", lenseq

		#If position belongs to the sequence select it
		if k in range(lenseq):
			tempSequences.append(sortedArray[s])
	tempSequences = np.array(tempSequences)
	#print "tempSequences", tempSequences.shape
	#print ""
	#Take the 10% of the sequences with the highest Initiation Rate
	topTen = int(len(tempSequences)*0.1)
	ninety = len(tempSequences) - topTen
	# print "topTen",topTen
	# print "ninety", ninety

	#Calculate the Pset
	As = Ts = Gs = Cs = 0
	for i in range(topTen):
		lenseq = len(tempSequences[i][0])

		if tempSequences[i][0][lenseq-(k+1)] == 'A':  #walk from right to left
			As += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'T':
			Ts += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'G':
			Gs += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'C':
			Cs += 1
		else:
			print "error top ten"

	pSetA = As/float(topTen)
	pSetT = Ts/float(topTen)
	pSetG = Gs/float(topTen)
	pSetC = Cs/float(topTen)
	# print "As", As
	# print "Ts", Ts
	# print "Gs", Gs
	# print "Cs", Cs
	# print "pSetA",pSetA
	# print "pSetT",pSetT
	# print "pSetG",pSetG
	# print "pSetC",pSetC
	
	#Calculate the Pback
	Ab = Tb = Gb = Cb = Nothingb = 0
	count = 0
	#print "range: ",topTen,len(tempSequences)
	#print "ninety", ninety
	for z in range(topTen,len(tempSequences)): #now check for the specific position in the rest 90%
		lenseq = len(tempSequences[z][0]) #length   !!!!PROSOXI!!!!!
		if tempSequences[z][0][lenseq-(k+1)] == 'A':  #walk from right to left
			Ab += 1
		elif tempSequences[z][0][lenseq-(k+1)] == 'T':
			Tb += 1
		elif tempSequences[z][0][lenseq-(k+1)] == 'G':
			Gb += 1
		elif tempSequences[z][0][lenseq-(k+1)] == 'C':
			Cb += 1
		else :
			print "error ninety"

	pBackA = Ab/float(ninety)
	pBackT = Tb/float(ninety)
	pBackG = Gb/float(ninety)
	pBackC = Cb/float(ninety)

	if int (pBackA+pBackT+pBackG+pBackC+0.0000000001) != 1:
		print "BIG ERROR", pBackA+pBackT+pBackG+pBackC

	
	#Now for each position calculate the Eset of each base
	if pBackA != 0 and (pSetA/pBackA) != 0:
		esetA = pSetA*math.log(pSetA/pBackA)
	else:
		esetA = 0
	arrayEsetA.append(esetA)
	#print "esetA", esetA


	if pBackT != 0 and (pSetT/pBackT) != 0:
		esetT = pSetT*math.log(pSetT/pBackT)
	else:
		esetT = 0
	arrayEsetT.append(esetT)


	if pBackG != 0 and (pSetG/pBackG) != 0:
		esetG = pSetG*math.log(pSetG/pBackG)
	else:
		esetG = 0
	arrayEsetG.append(esetG)

	if pBackC != 0 and (pSetC/pBackC) != 0:
		esetC = pSetC*math.log(pSetC/pBackC)
	else:
		esetC = 0
	arrayEsetC.append(esetC)

arrayEsetA = np.array(arrayEsetA)
arrayEsetT = np.array(arrayEsetT)
arrayEsetG = np.array(arrayEsetG)
arrayEsetC = np.array(arrayEsetC)
print "arrayEsetA", arrayEsetA.shape
print "arrayEsetT", arrayEsetT.shape
print "arrayEsetG", arrayEsetG.shape
print "arrayEsetC", arrayEsetC.shape
arrayEsetA = np.save('Ciandrini/arrayEsetA',arrayEsetA)
arrayEsetT = np.save('Ciandrini/arrayEsetT',arrayEsetT)
arrayEsetG = np.save('Ciandrini/arrayEsetG',arrayEsetG)
arrayEsetC = np.save('Ciandrini/arrayEsetC',arrayEsetC)

# # PLot
# barWidth = np.arange(len(arrayEsetA))
# plotWith4Subplots(4, arrayEsetA, 'A', arrayEsetT, 'T', arrayEsetG, 'G', arrayEsetC, 'C', barWidth, True, 'Conditional Entropy per base', 'Ciandrini/CEperBase.png')
# exit()
'''
'''
########################~~Summations~~########################################
#Now pass from every sequence and swap every base with its eset value
print "Start Summations"
# actualArrayWithNames = np.load('Gritsenko/forCE.npy')
# print "actualArrayWithNames loaded", actualArrayWithNames.shape
arrayEsetA = np.load('Ciandrini/arrayEsetA.npy')
arrayEsetT = np.load('Ciandrini/arrayEsetT.npy')
arrayEsetG = np.load('Ciandrini/arrayEsetG.npy')
arrayEsetC = np.load('Ciandrini/arrayEsetC.npy')

sumEsetSequences = []

for s in range(len(actualArrayWithNames)):
	lenseq = len(actualArrayWithNames[s][1]) # Length	!!!!PROSOXI!!!!!
	#print "lenseq", lenseq
	
	esetSequence = []
	for j in range(maxLength):

		#if the position is in between of the length of the sequence
		if j in range(lenseq):
			if actualArrayWithNames[s][1][j] == 'A':
				esetSequence.append(arrayEsetA[j])
			elif actualArrayWithNames[s][1][j] == 'T':
				esetSequence.append(arrayEsetT[j])
			elif actualArrayWithNames[s][1][j] == 'G':
				esetSequence.append(arrayEsetG[j])
			elif actualArrayWithNames[s][1][j] == 'C':
				esetSequence.append(arrayEsetC[j])
			else:
				print "error esetSequence" #esetSequence.append(arrayEsetNothing[j])

	Summation = sum(esetSequence)

	sumEsetSequences.append(Summation) #Append the Sum of the values of each Sequence-array.

sumEsetSequences = np.array(sumEsetSequences) # is an array that contains the Summation of the Eset values per sequence
print "Summation of the Conditional Entropy values is done: sumEsetSequences=", sumEsetSequences.shape

np.save("Ciandrini/sumEsetSequences", sumEsetSequences)
'''

sumEsetSequences = np.load("Ciandrini/sumEsetSequences.npy")
featureArray     = np.load("Ciandrini/featureArray.npy")
print "sumEsetSequences", sumEsetSequences.shape
print "featureArray",     featureArray.shape

tempArray = featureArray[:,0:-1]
tempArray = np.vstack((tempArray.T, sumEsetSequences)).T
featureArray = np.vstack((tempArray.T, featureArray[:,-1])).T
print "featureArray", featureArray.shape
print "featureArray", featureArray[0]
np.save('Ciandrini/featureArray', featureArray)
exit()

#print "sumEsetSequences", sumEsetSequences
# print "Done Conditional Entropy in ", time.time() - start_ConditionalEntropy_time

#plot sumEsetSequences against Init Rates

print type(sumEsetSequences[0])
for i in range(len(sumEsetSequences)):
	a = sumEsetSequences[i]
	print '{0:.10f}'.format(a)
exit()
y1 = sumEsetSequences
y1 ='{0:.10f}'.format(y1)
print y1
exit()
y2 = actualArrayWithNames[:,-1]
x  = np.arange(len(y1))
plotTwoScales(x, 'samples', y1, 'Summations', y2, 'InitRates', 'SummedValues from CE', 'Initiation Rates (log)', 'Gritsenko/plotSummations.png')


exit()
########

# Conditional Entropy per 3-mer

dimers,trimers,dimercounts,trimercounts = find_kmers('ATGC')

trimers = np.array(trimers) # 64 possibilities





essetTrimers = []

for k in range(maxLength-3):

	tempSequences = []

	#Select the sequences that have this length
	for s in range(len(sortedArray)):

		#Find Sequence's Length
		lenseq = int(sortedArray[s][2])
		#print "lenseq", lenseq

		#If position belongs to the sequence select it
		if k in range(lenseq):
			tempSequences.append(sortedArray[s])
	tempSequences = np.array(tempSequences)
	#print "tempSequences", tempSequences.shape
	#print ""
	#Take the 10% of the sequences with the highest Initiation Rate
	topTen = int(len(tempSequences)*0.1)
	ninety = len(tempSequences) - topTen
	# print "topTen",topTen
	# print "ninety", ninety


	#Calculate the Pset
	trisets = np.zeros(64)
	for i in range(topTen):
		lenseq = int(tempSequences[i][2])
		reverseseq = tempSequences[i][0][::-1]
		#print reverseseq
		#print lenseq, tempSequences[i][0]
		#print reverseseq[k:k+3]
		for tri_index in range(64):

			if reverseseq[k:k+3] == trimers[tri_index]:
				trisets[tri_index]+=1
	
	trisets = trisets / float(topTen)
	#print trisets
	
	tribacks = np.zeros(64)
	for i in range(topTen,len(tempSequences)):
		lenseq = int(tempSequences[i][2])
		reverseseq = tempSequences[i][0][::-1]
		#print reverseseq
		#print lenseq, tempSequences[i][0]
		#print reverseseq[k:k+3]
		for tri_index in range(64):

			if reverseseq[k:k+3] == trimers[tri_index]:
				tribacks[tri_index]+=1
	
	tribacks = tribacks / float(ninety)
	#print tribacks

	tri_essets = np.zeros(64)
	for p in range(64):
		if tribacks[p] != 0 and (trisets[p]/tribacks[p]) != 0:
			tri_essets[p] = trisets[p]*math.log(trisets[p]/tribacks[p])
		else:
			tri_essets[p] = 0
	#print tri_essets
	essetTrimers.append(tri_essets)

print len(essetTrimers)
print essetTrimers[0]

essetTrimers = np.array(essetTrimers)
#Plot the codon Conditional Entropies
for t_counter in range(0,64,4):
	print t_counter

	arrayEsset1 = essetTrimers[:,t_counter]
	arrayEsset2 = essetTrimers[:,t_counter+1]
	arrayEsset3 = essetTrimers[:,t_counter+2]
	arrayEsset4 = essetTrimers[:,t_counter+3]

	print "LENGTHS:",len(arrayEsset1),len(arrayEsset2),len(arrayEsset3),len(arrayEsset4)

	essetname1 = trimers[t_counter]
	essetname2 = trimers[t_counter+1]
	essetname3 = trimers[t_counter+2]
	essetname4 = trimers[t_counter+3]



	plotWith4Subplots(4, arrayEsset1, essetname1, arrayEsset2, essetname2, arrayEsset3, essetname3, arrayEsset4, essetname4, np.arange(len(arrayEsset1)), True, 'Conditional Entropy per Codon', 'Gritsenko/CEperAminoAcid'+`t_counter`+'.png')
