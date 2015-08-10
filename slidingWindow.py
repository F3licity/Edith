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
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
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

#Import My Own Functions
from myfunctions import parseGFFAttributes, parseGFF3, find_kmers, romanToNumeric, shortenSequence, findDuplcicates, searchIfInColumn, plotWith4Subplots, plotTwoScales, slidingWindow, findBasesFrequency, doRandmFrstRegression
######################################################################################################################################################

#Do the Slide Window Preparation
# featureArray = np.load('Gritsenko/featureArray.npy')
featureArray = np.load('Ciandrini/featureArray.npy')
print 'Loaded featureArray', featureArray.shape
# print featureArray[0]
# exit()
# print featureArray[:5,1]
# exit()
# print "", featureArray[:,3]
# exit()
# print featureArray[0,0]
# f = open('Ciandrini/mfe.txt','r')
# mfe = []
# for line in f:
#     mfe.append(line)
# mfe = np.array(mfe)
# mfe = mfe.astype('float')
# print "mfe", mfe.shape
# np.save('Ciandrini/mfe', mfe)
# exit()
# for i in range(len(mfe)):
# 	featureArray[i][89] = mfe[i]
# np.save('Gritsenko/featureArray.npy', featureArray)

#### Sliding Window - Prepare Features from the Fragments - forRegression array ####
windowSize = 9
step = 1
fragmentedSeqs = []
for i in range(len(featureArray)):

	# reverseseq = featureArray[i][0][::-1]   # cause I want all the sequences alligned on the 3'
	reverseseq = featureArray[i][1][::-1]
	chunks = slidingWindow(reverseseq, windowSize, step)
	chunks = list(chunks)
	fragmentedSeqs.append(chunks)
#Do not save this array, you cant load it back more than a window.


cs=[]
for i in range(len(featureArray)):
	#print "Sequence of length", featureArray[i][2]
	c = 0
	for j in fragmentedSeqs[i]:
		
		c +=1  #Count Fragments
	#print c
	cs.append(c)
cs = np.array(cs)
maxGroups=cs.max(axis=0)
minGroups=cs.min(axis=0)
print "maxGroups", maxGroups
print "minGroups", minGroups
maximum = (4 + (maxGroups*84))
print "maximum", maximum

# upstreamsequences = np.load('Gritsenko/upstreamsequences.npy')

# Max No of fragments is 19 coming from the length boundary 140bp

allfeaturenames = []
basenames = ["A","T","G","C"]

forRegression = []
for i in range(len(featureArray)):
	#For Gritsenko
	length   = featureArray[i][3]
	stress   = featureArray[i][87]
	ce       = featureArray[i][88]
	mfe      = featureArray[i][89]
	initRate = featureArray[i][95]

	# #For Ciandrini
	# length     = featureArray[i][3]
	# stress     = featureArray[i][88]
	# ce         = featureArray[i][90]
	# mfe        = featureArray[i][89]
	# initRate   = featureArray[i][91]


	if (len(fragmentedSeqs[i]) == maxGroups):
		allfeaturenames = ["Length","Stress","CE","MFE"]

	#print "Length", featureArray[i][2], existingGenes[i][2], upstreamsequences[i][1]

	#Get Features of the whole sequence
	line = np.hstack((length, stress, ce, mfe))

	# print len(featureArray[i][0])
	# #Get Features of the Fragments
	
	
	
	count = 0
	for j in fragmentedSeqs[i]:
		temp = []
		something = []
		# j is the fragment
		count +=1
		# print j
		
		
		#Bases Frequencies
		#findBasesFrequency(j) = Afrequency, Tfrequency, Gfrequency, Cfrequency
		temp = findBasesFrequency(j)
		#line = np.hstack((line, temp))

		# print "temp", temp
		# print "temp size", len(temp)
		# print ""
		dimers,trimers,dimercounts,trimercounts = find_kmers(j)

		# print "dimers",dimers
		# print "j",j
		# print "dimercounts",dimercounts
		# exit()

		maxNoDimers = float(windowSize-1)   # The length of each window is 10
		maxNoTrimers = float(windowSize-2)

		for dim in dimercounts:
			if maxNoDimers != 0:
				something.append(dim/maxNoDimers)
			else:
				something.append(dim)
		for trim in trimercounts:
			if maxNoTrimers != 0:
				something.append(trim/maxNoTrimers)
			else:
				something.append(trim)
		# print "something", len(something)
		# print ""
		line = np.hstack((line, temp, something))
		if (len(fragmentedSeqs[i]) == maxGroups):
			basenames2 = [s + " "+`count` for s in basenames]
			dimers2 = [s + " "+`count` for s in dimers]
			trimers2 = [s + " "+`count` for s in trimers]
			allfeaturenames = np.hstack((allfeaturenames,basenames2,dimers2,trimers2))
		# print "temp, something", temp, something
		
		# print "line", line
	# print count
	# print "line **", len(line)

	# # Max no of fragments is 19. So a line must have 4+(19x84)=1600 (4+16+64=84)
	# difference = 1600 - (4+(count*84))
	# Max no of fragments is 19. So for only the bases frequency mus have 4+(19*4)=80
	difference = maximum - (4+(count*84))
	lst = [-1] * difference

	line = np.hstack((line, lst))
	line = np.hstack((line,initRate))
	# print len(line)
	forRegression.append(line)

forRegression = np.array(forRegression)
print "forRegression", forRegression.shape
# print "forRegression", forRegression[0]
np.save('Ciandrini/forRegression', forRegression)

forRegression = np.load('Ciandrini/forRegression.npy')
print "forRegression", forRegression.shape


X = forRegression[:,0:-1].astype(float) #Features
y = forRegression[:,-1].astype(float) #Target

print "len(X)", len(X)
print "len(X.T)", len(X.T)


pearsonsCorrelations = []
spearmanCorrelations = []
for i in range(len(X.T)): #per column.

	indexes = np.where(X[:,i] > -1)
	indexes = indexes[0]

	pC = scipy.stats.pearsonr(X[indexes,i], y[indexes])
	pearsonsCorrelations.append(pC)

	sC = scipy.stats.spearmanr(X[indexes,i], y[indexes])
	spearmanCorrelations.append(sC)


pearsonsCorrelations = np.array(pearsonsCorrelations)
spearmanCorrelations = np.array(spearmanCorrelations)
print "pearson", pearsonsCorrelations.shape
print "spearman",spearmanCorrelations.shape


h = []
p = []
z = []
q = []
x = []

for i in range(len(spearmanCorrelations)):
	if math.fabs(spearmanCorrelations[i][0]) > 0.01:
		h.append(spearmanCorrelations[i][0]) # Spearman correlation coefficient
		p.append(spearmanCorrelations[i][1]) # Spearman P values
		x.append(i)
	# if math.fabs(pearsonsCorrelations[i][0]) > 0.03:
	# 	z.append(pearsonsCorrelations[i][0]) # Pearson Correlation
	# 	q.append(pearsonsCorrelations[i][1]) # Pearson P values

# ind = np.arange(len(h)) #width of a bar

# host = host_subplot(111, axes_class=AA.Axes)
# par1 = host.twinx()

# host.set_xlabel('Features')
# host.set_ylabel('Pearson Cor. Coefficient')
# par1.set_ylabel('P-Value')

# p1 = host.bar(x, h, color='g', alpha=0.5, linewidth=0)
# p2 = par1.bar(x, p, color='y', alpha=0.5, linewidth=0)

# host.set_xticks(x)
# host.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
# host.axis["left"].label.set_color('green')
# par1.axis["right"].label.set_color('yellow')
# plt.draw()
# plt.show()


#plt.savefig(filename)



plt.clf()
ind = np.arange(len(x)) #width of a bar

f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
#f.tight_layout()
ax1.bar(ind, h, color='g', alpha=0.9, linewidth=0)
# ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
ax1.set_title('Spearman Cor.Coefficient')
ax1.set_xticks(ind+0.4)
ax1.set_xticklabels(x)
# ax1.grid()
ax2.bar(ind, p, color='r', alpha=0.9, linewidth=0)
# ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
ax2.set_title('P value')
ax2.set_xticks(ind+0.4)
ax2.set_xticklabels(allfeaturenames[x])
# ax2.grid()
f.suptitle('Feats with higher Correlation than 0.01')
f.savefig('Ciandrini/StupidFeatsSelectionSpearman.png')

exit()
# plt.show()



# #DO RANDOM FOREST REgression
X = forRegression[:,x].astype(float) #Features
y = forRegression[:,-1].astype(float) #Target

skf = cross_validation.KFold(len(y),n_folds=5)
# skf = cross_validation.StratifiedKFold(b,n_folds=5)

# saveImageName ='Gritsenko/SlidingWindow/rndmForest'
# correlationImageName = 'Gritsenko/SlidingWindow/RndmForestCorrelationFold'

# doRandmFrstRegression(X, y, skf, saveImageName, correlationImageName)

saveImageName ='Ciandrini/GradientBoosting/GBReg'
correlationImageName = 'Ciandrini/GradientBoosting/GBReg'

#clf = linear_model.Lasso(alpha=0.1, normalize=True)

doRandmFrstRegression(X, y, skf, saveImageName, correlationImageName)