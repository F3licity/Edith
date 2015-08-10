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
arrayNames        = np.load("arrayNames.npy")
arrayFeatures     = np.load("arrayFeatures.npy")
codingsLengths     = np.load("codingsLengths.npy")
arrayCodSequences = np.load("arrayCodSequences.npy")

print "arrayNames : ", arrayNames.shape
print "arrayFeatures : ", arrayFeatures.shape
print "codingsLengths : ", codingsLengths.shape
print "arrayCodSequences : ",arrayCodSequences.shape

X = codingsLengths
y = arrayFeatures[:,-1]

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
    X_test = X_test.reshape( (len(X_test), 1) )
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
