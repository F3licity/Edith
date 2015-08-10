###########################################################
#Ordinary Least Squares

actualArray = np.array(actualArray)
print "shape",actualArray.shape

clf = linear_model.LinearRegression()
clf.fit ([actualArray[0], actualArray[1]], [0, 1])
print "Least Squares coefficient", clf.coef_

#Define x and y
X = actualArray[:,:-1] #np.array([actualArray[:,0:4]]).T
print X.T.shape
y = actualArray[:,-1]
print y.shape
# print X 
# print y
X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, y, test_size=0.2)

# print "x train is ", X_train
# print "y train is ", y_train

# print "x test is ", X_test
# print "y test is ", y_test

model = linear_model.LinearRegression()
#print model
model.fit(X_train, y_train)
#print model

predicted = model.predict(X_test)

plt.plot(predicted, 'r')
plt.plot(y_test,'g')
plt.title("LinearRegression Only Length")
plt.xlabel("features axis")
plt.ylabel("Init. rates axis")
#plt.savefig('LinearRegression Only Length.png')
#plt.show()
plt.clf()
print('MeanAbsoluteError: {}'.format(metrics.mean_absolute_error(model.predict(X_test), y_test)))
print "--------------------------------------------------------------------------------------------------end of Least Squares"

# ----{Random Forest}----
# class sklearn.ensemble.RandomForestRegressor(n_estimators=10, criterion='mse', max_depth=None,
#     min_samples_split=2, min_samples_leaf=1, max_features='auto', max_leaf_nodes=None, bootstrap=True,
#     oob_score=False, n_jobs=1, random_state=None, verbose=0, min_density=None, compute_importances=None


res = RandomForestRegressor(n_estimators=100)
#res = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
res = res.fit(X_train, y_train)

Y_predicted = res.predict(X)
Y_test_predicted = res.predict(X_test)

print "params", res.get_params([True])
print "X is ", X
print "Prediction of X is ", Y_predicted

print "Random Forest score", res.score(X_test, y_test)

print "MSE:", mean_squared_error(y_test, Y_test_predicted)
# plt.plot(res.predict(X_test),'r')
# plt.plot(y_test,'g')
# plt.title("RandomForest Only Length")

# plt.xlabel("Length")
# plt.ylabel(" Error")
# plt.title(" error of predictions")
# plt.scatter(X_test[:,0], y_test-Y_test_predicted , color='b', alpha=0.5)
# plt.savefig("absolute_error_test.png")
# plt.clf()

# plt.xlabel("length axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Length")
# plt.scatter(actualArray[:,0], y, color='g', alpha=0.5)
# plt.scatter(actualArray[:,0], Y_predicted, color='r', alpha=0.5)
# plt.savefig("Scatter against Length_test.png")
# plt.clf()

# plt.xlabel("A axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against As")
# plt.scatter(actualArray[:,1], y, alpha=0.5)
# plt.savefig("Scatter against As_test.png")
# plt.clf()

# plt.xlabel("T axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Ts")
# plt.scatter(actualArray[:,2], y, alpha=0.5)
# plt.savefig("Scatter against Ts_test.png")
# plt.clf()

# plt.xlabel("G axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Gs")
# plt.scatter(actualArray[:,3], y, alpha=0.5)
# plt.savefig("Scatter against Gs_test.png")
# plt.clf()

# plt.xlabel("C axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Cs")
# plt.scatter(actualArray[:,4], y, alpha=0.5)
# plt.savefig("Scatter against Cs_test.png")
# plt.clf()


#plt.show()
#plt.clf()
print "--------------------------------------------------------------------------------------------------end of RandomForest"


        #empties = ""
        #for x in range(maxLength-lenseq):
        #    empties = empties+"_"
        #empties = ["" for x in range(maxLength-lenseq)]
        #print "empties", empties
        #actualArray[i][0] = actualArray[i][0] + empties
        #np.concatenate((actualArray[i][0],empties), axis=0)
        #print actualArray[i][0]
----------------------------------------------------------------------------------------------------------------------------------

#Heatmap
data = Data([
    Heatmap(
        z=actualArray_scaled[:,0] # 1 length
    )
])
plot_url = py.plot(data, filename='Length' + `csvFile`+ ' heatmap')
data = Data([
    Heatmap(
        z=actualArray_scaled[:,85]
    )
])
plot_url = py.plot(data, filename='CE summations ' + `csvFile`+ ' heatmap')

PLots

# ind = np.arange(len(arrayEsetA)) #is about the width of a bar
# ####subplot(nrows, ncols, plot_number)###
# plt.close('all')
# plt.figure()
# # plt.subplots_adjust(hspace=.5, wspace=.5)
# plt.subplots(nrows=5, ncols=1, sharex=True, sharey=True)
# plt.tight_layout()
# plt.subplot(5,1,1)
# plt.bar(ind, arrayEsetA,0.8,color='b',linewidth=0)
# plt.title('Conditional Entropy A')

# plt.subplot(5,1,2) # second subplot in the figure
# plt.bar(ind, arrayEsetT,0.8,color='r',linewidth=0)
# plt.title('Conditional Entropy T')

# plt.subplot(5,1,3)
# plt.bar(ind, arrayEsetG,0.8,color='g',linewidth=0)
# plt.title('Conditional Entropy G')

# plt.subplot(5,1,4)
# plt.bar(ind, arrayEsetC,0.8,color='y',linewidth=0)
# plt.title('Conditional Entropy C')

# plt.subplot(5,1,5)
# plt.bar(ind, arrayEsetNothing,0.8,color='k',linewidth=0)
# plt.title('Conditional Entropy No-base')

# plt.xlabel('Direction: 3 to 5')
# plt.suptitle('Conditional Entropy per base for ' + `csvFile`)
# # plt.savefig('Conditional Entropy per base for ' + `csvFile`+'.png')
# plt.show()
# plt.clf()
___________________________________________________________________________
Latex file
\begin{figure}[!ht]
\begin{center}
\includegraphics[scale=0.8]{bayes.eps}
\end{center}
\caption{Bayesian network with use of Bioinformatics data sources. One 
can see that the nodes may represent different variables coming from 
different sources. \label{bayes}}
\end{figure}
_____________________________________________________________________________

# #Feature Selection
# slctKbest = SelectKBest(score_func = f_regression, k = 10)
# slctKbest.fit_transform(X,y)    # not just fit, in order to get the kBest
# print "get_params", slctKbest.get_params(deep = True)    #If True, will return the parameters for this estimator and contained subobjects that are estimators.
# print "getSupport", slctKbest.get_support(indices = True)
# # get_params {'k': 10, 'score_func': <function f_regression at 0x39e2320>}
# # getSupport [ 4  9 16 26 33 38 41 43 64 80]    Those are my informative features
# print "slctKbest.scores_", slctKbest.scores_
# print "slctKbest.pvalues_", slctKbest.pvalues_

# exit()
# newFeatures = slctKbest.fit_transform(X,y)


###################################################################################################################################
for degree in [1]:
    svr_linear = svm.SVR(kernel='linear', C=1, cache_size=200) #SupportVectorMachine - Support Vector Regression
    svr_rbf = svm.SVR(kernel='rbf', C=1, gamma=0.0, degree=degree, cache_size=200)
    svr_poly = svm.SVR(kernel='poly', C=1, degree=degree, coef0=0.0, cache_size=200)
    svr_sigmoid = svm.SVR(kernel='sigmoid', C=1, degree=degree, coef0=0.0, cache_size=200)

#class sklearn.svm.SVR(kernel='rbf', degree=3, gamma=0.0, coef0=0.0, tol=0.001, C=1.0,
# epsilon=0.1, shrinking=True, probability=False, cache_size=200, verbose=False, 
#max_iter=-1, random_state=None)

    #Get the score for each Fold 
    #SupportVectorRegression(Model, data, Target, Folds)
    scores_linear = cross_validation.cross_val_score(svr_linear, X, y, cv=5)
    print "Cross Validation Scores Linear", scores_linear
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores_linear.mean(), scores_linear.std()))

    scores_rbf = cross_validation.cross_val_score(svr_rbf, X, y, cv=5)
    print "Cross Validation Scores RBF", scores_rbf
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores_rbf.mean(), scores_rbf.std()))

    scores_poly = cross_validation.cross_val_score(svr_poly, X, y, cv=5)
    print "Cross Validation Scores Polynomial", scores_poly
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores_poly.mean(), scores_poly.std()))

    scores_sigmoid = cross_validation.cross_val_score(svr_sigmoid, X, y, cv=5)
    print "Cross Validation Scores sigmoid", scores_sigmoid
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores_sigmoid.mean(), scores_sigmoid.std()))

y_linear = svr_linear.fit(X, y).predict(X)
y_rbf = svr_rbf.fit(X, y).predict(X)
y_poly = svr_poly.fit(X, y).predict(X)
y_sigmoid = svr_sigmoid.fit(X, y).predict(X)

pr=np.arange(len(X[:,1]))
plt.scatter(pr, y, c='k', label='data')
plt.hold('on')
plt.plot(pr, y_linear, c='r', label='Linear model')
plt.plot(pr, y_rbf, c='g', label='RBF model')
plt.plot(pr, y_poly, c='b', label='Polynomial model')
plt.plot(pr,y_sigmoid, c='orange', label='sigmoid')
plt.xlabel('#of rows')
plt.ylabel('initiation rates (log)')
plt.title('Support Vector Regression - Length')
plt.legend()
plt.savefig("Support Vector Regression - Length.png")
# plt.show()
plt.clf()

absoluteDifference = np.absolute(y-y_linear)
pr=np.arange(len(X[:,1]))
plt.scatter(pr, absoluteDifference, c='k', label='data')
plt.savefig("Absolute Error - Linear SVM.png")
# plt.show()
plt.clf()


plt.scatter(X[:,1], y, c='k', label='data')
plt.hold('on')
plt.plot(pr, y_linear, c='r', label='Linear model')
# plt.plot(X[:,1], y_rbf, c='g', label='RBF model')
# plt.plot(X[:,1], y_poly, c='b', label='Polynomial model')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression - A')
plt.legend()
plt.savefig("Support Vector Regression - A.png")
#plt.show()
plt.clf()

plt.scatter(X[:,2], y, c='k', label='data')
plt.hold('on')
plt.plot(X[:,2], y_linear, c='r', label='Linear model')
# plt.plot(X[:,2], y_rbf, c='g', label='RBF model')
# plt.plot(X[:,2], y_poly, c='b', label='Polynomial model')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression - T')
plt.legend()
plt.savefig("Support Vector Regression - T.png")
#plt.show()
plt.clf()

plt.scatter(X[:,3], y, c='k', label='data')
plt.hold('on')
plt.plot(X[:,3], y_linear, c='r', label='Linear model')
# plt.plot(X[:,3], y_rbf, c='g', label='RBF model')
# plt.plot(X[:,3], y_poly, c='b', label='Polynomial model')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression - G')
plt.legend()
plt.savefig("Support Vector Regression - G.png")
#plt.show()
plt.clf()

plt.scatter(X[:,4], y, c='k', label='data')
plt.hold('on')
plt.plot(X[:,4], y_linear, c='r', label='Linear model')
# plt.plot(X[:,4], y_rbf, c='g', label='RBF model')
# plt.plot(X[:,4], y_poly, c='b', label='Polynomial model')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression - C')
plt.legend()
plt.savefig("Support Vector Regression - C.png")
#plt.show()
plt.clf()

# plt.xlabel("Length")
# plt.ylabel(" Error")
# plt.title(" error of predictions")
# plt.scatter(X_test[:,0], y_test-Y_test_predicted , color='b', alpha=0.5)
# plt.savefig("absolute_error_test.png")
# plt.clf()

#clf.fit(X, y) 
# SVR(C=1, cache_size=200, coef0=0.0, degree=3, epsilon=0.1, gamma=0.0,
#   kernel='linear', max_iter=-1, probability=False, random_state=None,
#   shrinking=True, tol=0.001, verbose=False)
#------------------------------------------------------------------------------
#  R^2 = (1-u/v) where u =  ((y_true - y_pred) ** 2).sum()
# and v =  ((y_true - y_true.mean()) ** 2).sum()
# best value of R^2 =1 
# r_squared = clf.score(X, y, sample_weight=None)
# print "R-Squared", r_squared
# I don't use this one because it calculateds the Squared error relying on all the data, training & test (bias)
#-------------------------------------------------------------------------------





# plt.xlabel("length axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Length")
# plt.scatter(actualArray[:,0], y, color='g', alpha=0.5)
# plt.scatter(actualArray[:,0], Y_predicted, color='r', alpha=0.5)
# plt.savefig("Scatter against Length_test.png")
# plt.clf()

# plt.xlabel("A axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against As")
# plt.scatter(actualArray[:,1], y, alpha=0.5)
# plt.savefig("Scatter against As_test.png")
# plt.clf()

# plt.xlabel("T axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Ts")
# plt.scatter(actualArray[:,2], y, alpha=0.5)
# plt.savefig("Scatter against Ts_test.png")
# plt.clf()

# plt.xlabel("G axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Gs")
# plt.scatter(actualArray[:,3], y, alpha=0.5)
# plt.savefig("Scatter against Gs_test.png")
# plt.clf()

# plt.xlabel("C axis")
# plt.ylabel("Init. rates axis")
# plt.title("Scatter against Cs")
# plt.scatter(actualArray[:,4], y, alpha=0.5)
# plt.savefig("Scatter against Cs_test.png")
# plt.clf()
print "--------------------------------------------------------------------------------------------------end of SVM, R-Squared"


#Make a plot
plt.clf()
#~~Make a plot of Initiation Rates~~
plt.close('all')
plt.figure()
# plt.subplot(1,2,1)
y = sortedArray[:,-1] #Initiation rates
x = np.arange(len(y))
plt.plot(x,y,'g', label='Alexey')
plt.xlabel("Indices")
plt.ylabel("Init Rates (log)")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.grid(True)
plt.savefig("InitRatesAlexeyThreshold.png")
exit()#Correlation Fitness, Segments, AvgFitness
pearsonsCorrelations = []
spearmanCorrelations = []
feats = featureArray[:,-4:-1]
Target = featureArray[:,-1]
# SUbplots  COrrelation - P value
for i in range(len(feats.T)): #per column.
    pC = scipy.stats.pearsonr(feats[:,i], Target)
    pearsonsCorrelations.append(pC)
    # print pC
    sC = scipy.stats.spearmanr(feats[:,i], Target)
    spearmanCorrelations.append(sC)
    # print sC

pearsonsCorrelations = np.array(pearsonsCorrelations)
spearmanCorrelations = np.array(spearmanCorrelations)
# print "pearsonsCorrelations", pearsonsCorrelations.shape
# print "spearmanCorrelations", spearmanCorrelations.shape

h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
z = spearmanCorrelations[:,0] # Spearman
x = np.arange(len(feats.T)) # Features
p = pearsonsCorrelations[:,1] # Pearson P values
q = spearmanCorrelations[:,1] # Spearman P values

plt.clf()
ind = np.arange(len(feats.T)) #width of a bar
f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
f.tight_layout()
ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
ax1.set_xticks(np.arange(len(feats.T)))
ax1.grid()
ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
ax2.set_title('P value')
ax2.set_xticks(np.arange(len(feats.T)))
ax2.grid()
f.suptitle('Correlation avg Fitness')
f.savefig('Corr-Pvalue-AvgFitness.svg')
#################################################################################################
15-07-2015.
#Correlation Fitness, Segments, AvgFitness
pearsonsCorrelations = []
spearmanCorrelations = []
feats = featureArray[:,-4:-1]
Target = featureArray[:,-1]
# SUbplots  COrrelation - P value
for i in range(len(feats.T)): #per column.
    pC = scipy.stats.pearsonr(feats[:,i], Target)
    pearsonsCorrelations.append(pC)
    # print pC
    sC = scipy.stats.spearmanr(feats[:,i], Target)
    spearmanCorrelations.append(sC)
    # print sC

pearsonsCorrelations = np.array(pearsonsCorrelations)
spearmanCorrelations = np.array(spearmanCorrelations)
# print "pearsonsCorrelations", pearsonsCorrelations.shape
# print "spearmanCorrelations", spearmanCorrelations.shape

h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
z = spearmanCorrelations[:,0] # Spearman
x = np.arange(len(feats.T)) # Features
p = pearsonsCorrelations[:,1] # Pearson P values
q = spearmanCorrelations[:,1] # Spearman P values

plt.clf()
ind = np.arange(len(feats.T)) #width of a bar
f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
f.tight_layout()
ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
ax1.set_xticks(np.arange(len(feats.T)))
ax1.grid()
ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
ax2.set_title('P value')
ax2.set_xticks(np.arange(len(feats.T)))
ax2.grid()
f.suptitle('Correlation avg Fitness')
f.savefig('Corr-Pvalue-AvgFitness.svg')
-------
plt.clf()
#Plot the Fitness against the Init. Rates
x1 = np.log( sortedFeaturesArray[:,-5] )# mRNA counts
print "mRNA", x1
x2 = np.log(sortedFeaturesArray[:,-4] )# Ribosome counts
print "Ribosome", x2

x3 = sortedFeaturesArray[:,-3] # Fitness      Not anymore in this column...
x4 = sortedFeaturesArray[:,-2] # No of Segments
y  = sortedFeaturesArray[:,-1] # Initiation Rates
x  = np.arange(len(y))
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2) #, sharex='col', sharey='row')
ax1.scatter(x, y, c='g',s=0.4,lw = 0)
ax1.scatter(x,x1,c='r',s=1.5,lw = 0)
ax1.set_title('mRNA counts')
ax2.scatter(x, y, c='g',s=0.4,lw = 0)
ax2.scatter(x,x2,c='r',s=1.5,lw = 0)
ax2.set_title('Ribosome counts')
ax3.scatter(x, y, c='g',s=0.4,lw = 0)
ax3.scatter(x,x3,c='r',s=1.5,lw = 0)
ax3.set_title('Fitness')
ax4.scatter(x, y, c='g',s=0.4,lw = 0)
ax4.scatter(x,x4,c='r',s=1.5,lw = 0)
ax4.set_title('NoSegments')
plt.savefig('Stress/ALLFEats/4plotsExtraNumbers.svg')
# plt.show()
# y = []
# for i in range(len(sortedFeaturesArray)):
#     y.append(sortedFeaturesArray[i][-3]/sortedFeaturesArray[i][-2])
# x  = np.arange(len(y))
# plt.figure()
# plt.scatter(x,y,c='r',s=1.5,lw = 0, label='Fitness/NoofSegments')
# plt.scatter(x, sortedFeaturesArray[:,-1], c='g',s=1,lw = 0, label='InitRates')
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=2, mode="expand", borderaxespad=0.)
# plt.grid(True)
# # plt.show()
# plt.savefig('Stress/ALLFEats/PlotofAvgFitness.svg')

USEFULLL
#Correlation
pearsonsCorrelations = []
spearmanCorrelations = []
rawFeatures = actualArrayWithNames[:,2:-1]
rawTarget = actualArrayWithNames[:,-1]
Features = actualArrayWithNames[:,2:-1].astype(float)
Target = actualArrayWithNames[:,-1].astype(float)
print "Features", Features.shape
'''
# SUbplots  COrrelation - P value
# for i in range(len(Features.T)): #per column.
#     pC = scipy.stats.pearsonr(Features[:,i], Target)
#     pearsonsCorrelations.append(pC)
#     # print pC
#     sC = scipy.stats.spearmanr(Features[:,i], Target)
#     spearmanCorrelations.append(sC)
#     # print sC

# pearsonsCorrelations = np.array(pearsonsCorrelations)
# spearmanCorrelations = np.array(spearmanCorrelations)
# # print "pearsonsCorrelations", pearsonsCorrelations.shape
# # print "spearmanCorrelations", spearmanCorrelations.shape

# h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
# z = spearmanCorrelations[:,0] # Spearman
# x = np.arange(len(Features.T)) # Features
# p = pearsonsCorrelations[:,1] # Pearson P values
# q = spearmanCorrelations[:,1] # Spearman P values

# plt.clf()
# ind = np.arange(len(Features.T)) #width of a bar
# f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
# f.tight_layout()
# ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
# ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
# ax1.set_xticks(np.arange(len(Features.T)))
# ax1.grid()
# ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
# ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
# ax2.set_title('P value')
# ax2.set_xticks(np.arange(len(Features.T)))
# ax2.grid()
# f.suptitle('Feature Correlation at fold')
# f.savefig('Corr-Pvalue-Test.svg')


#PLots of Correlation and P value SeparateImages
for i in range(len(Features.T)): #per column.
    pC = scipy.stats.pearsonr(Features[:,i], Target)
    pearsonsCorrelations.append(pC)
    # print pC
    sC = scipy.stats.spearmanr(Features[:,i], Target)
    spearmanCorrelations.append(sC)
    # print sC

pearsonsCorrelations = np.array(pearsonsCorrelations)
spearmanCorrelations = np.array(spearmanCorrelations)
# print "pearsonsCorrelations", pearsonsCorrelations.shape
# print "spearmanCorrelations", spearmanCorrelations.shape

h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
z = spearmanCorrelations[:,0] # Spearman
x = np.arange(len(Features.T)) # Features

plt.clf()
ind = np.arange(len(Features.T))
plt.figure(figsize=(26, 12))
pc = plt.bar(ind, h, color='g', alpha=0.5)
sc = plt.bar(ind, z, color='r', alpha=0.5)

plt.title("Correlation Stress AgGregation")
# plt.bar(range(len(stressedFeatureValues.T)), h,
#        color="g", align="center")
plt.xticks(range(len(Features.T)))
plt.xlim([0, len(Features.T)])
plt.grid(True)
art=[]
plt.legend( (pc, sc), ('Pearson', 'Spearman') )
#plt.show()
plt.savefig("Stress/ALLStress/CorrelationStress"+ `csvFile` + ".png", additional_artists=art, bbox_inches='tight')
plt.clf()


h = pearsonsCorrelations[:,1] # Pearson P values
z = spearmanCorrelations[:,1] # Spearman
x = np.arange(len(Features.T)) # Features

ind = np.arange(len(Features.T))
plt.figure(figsize=(26, 12))
pc = plt.bar(ind, h, color='g', alpha=0.5)
sc = plt.bar(ind, z, color='r', alpha=0.5)

plt.title("P Values Stress Aggregation")
# plt.bar(range(len(stressedFeatureValues.T)), h,
#        color="g", align="center")
plt.xticks(range(len(Features.T)))
plt.xlim([0, len(Features.T)])
plt.grid(True)
art=[]
plt.legend( (pc, sc), ('Pearson', 'Spearman') )
#plt.show()
plt.savefig("Stress/ALLStress/PValues"+ `csvFile` + ".png", additional_artists=art, bbox_inches='tight')
plt.clf()
'''

########################################
#Stress
# for j in range(len(stressedGenes)):
#     for i in range(noRows): # no or Rows in ActualArray = Existing Genes


#         if (actualArrayWithNames[i][0] == stressedGenes[j][0]) and (actualArrayWithNames[i][-2] == '-2'):
#             actualArrayWithNames[i][-2] = stressedGenes[j][2]
#         elif (actualArrayWithNames[i][0] == stressedGenes[j][0]) and (actualArrayWithNames[i][-2] != '-2'):
#             extraLine = deepcopy(actualArrayWithNames[i])
            
#             extraLine[-2] = stressedGenes[j][2]
#             actualArrayWithNames = np.vstack((actualArrayWithNames,extraLine))

# print "actualArrayWithNames", actualArrayWithNames.shape

# stressedFeatureValues = np.vstack((a.T,b)).T
# print "stressedFeatureValues", stressedFeatureValues.shape


# ~~~~PLOt  the Stress -----
# plt.clf()

# y = sortedArray[:,-2].astype(int)#Stress
# # print "Stress", y
# # maxNo = y.max(axis=0)
# # minNo = y.min(axis=0)
# # print "maxNo", maxNo
# # print "minNo", minNo

# z = sortedArray[:,-1].astype(float) #Initiation Rates


# x = np.arange(len(y))
# plt.figure(figsize=(26, 12))
# plt.plot(y,z,'-o')
# plt.savefig("Stress/ALLFEats/PlotDot Stress-InitRates for " + `csvFile` + ".svg")
# exit()

# plt.plot(x,y,'g', label='Aggregation of Stress Terms')
# plt.plot(x,z,'r', label='InitiationRates')
# plt.xlabel("Indices")
# plt.ylabel("Stress")
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=2, mode="expand", borderaxespad=0.)
# # plt.show()
# #plt.savefig("Stress/ALLFEats/Plot Stress-InitRates for " + `csvFile` + ".png")
# plt.clf()
# exit()

######################################################################################################################################
#####################################################################################################################################

'''
# Do REgression So Far 86 Features
X = rawFeatures #a #stressedFeatureValues  # I should use Features and not rawFeatures
y = rawTarget
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
    print("TRAIN:", train_index, "TEST:", test_index)
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
    plt.suptitle('rndmForest_stress'+`fold`+'Fold')
    plt.savefig('Stress/ALLStress/rndmForest'+ `csvFile` +'Fold '+`fold`+'.png', additional_artists=art, bbox_inches='tight')
    # plt.show()
    plt.clf()

    #Calculate the Correlation at eeach foldach fold
    pearsonsCorrelations = []
    spearmanCorrelations = []
    
    for i in range(len(rawFeatures.T)): #per column.   # I should use Features and not rawFeatures
        pC = scipy.stats.pearsonr(rawFeatures[:,i], Target)
        pearsonsCorrelations.append(pC)
        # prawrint pC
        sC = scipy.stats.spearmanr(rawFeatures[:,i], Target)
        spearmanCorrelations.append(sC)
        # print sC

    pearsonsCorrelations = np.array(pearsonsCorrelations)
    spearmanCorrelations = np.array(spearmanCorrelations)
    # print "pearsonsCorrelations", pearsonsCorrelations.shape
    # print "spearmanCorrelations", spearmanCorrelations.shape

    h = pearsonsCorrelations[:,0] # Pearson correlation coefficient
    z = spearmanCorrelations[:,0] # Spearman
    x = np.arange(len(rawFeatures.T)) # Features
    p = pearsonsCorrelations[:,1] # Pearson P values
    q = spearmanCorrelations[:,1] # Spearman P values

    plt.clf()
    ind = np.arange(len(rawFeatures.T)) #width of a bar
    f, (ax1, ax2) = plt.subplots(2, sharex = True, figsize=(25,15))
    f.tight_layout()
    ax1.bar(ind, h, color='g', alpha=0.5, linewidth=0)
    ax1.bar(ind, z, color='r', alpha=0.5, linewidth=0)
    ax1.set_xticks(np.arange(len(rawFeatures.T)))
    ax1.grid()
    ax2.bar(ind, p, color='g', alpha=0.5, linewidth=0)
    ax2.bar(ind, q, color='r', alpha=0.5, linewidth=0)
    ax2.set_title('P value')
    ax2.set_xticks(np.arange(len(rawFeatures.T)))
    ax2.grid()
    f.suptitle('Feature Correlation at fold '+ `fold`)
    f.savefig('Stress/ALLStress/RndmForestCorrelationFold'+ `fold`+'.svg')
    plt.clf()
    fold += 1

# print "Prediction of X is ", Y_predicted
allScores = np.array(allScores)
print "Random Forest Scores:", allScores
print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(y_predictTrain, y_train)))
print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(y_predictTest, y_test)))
exit()
'''

##############################################################################################
#############~~Position Depending~~ ################################
start_ConditionalEntropy_time = time.time()
print "start conditional Entropy"
# np.save('Gritsenko/sortedArray', sortedArray)
sortedArray = np.load('Gritsenko/sortedArray.npy')
print "sortedArray  loaded, is the actualArrayWithNames Array (without the names) sorted by Init Rates: ", sortedArray.shape
print ""

temp = sortedArray[:,2].astype(int) # LEnGTH
maxLength=temp.max(axis=0)
print "maxLength", maxLength # Find max length
minLength=temp.min(axis=0)
print "minLength", minLength
print ""

#Calculate the Conditional Entropy (=Eset)
#Take the 10% of the sequences with the highest Initiation Rate
topTen = int(len(sortedArray)*0.1)
ninety = len(sortedArray) - topTen
print "90%% of the array is %d sequences" % ninety

#sortedArray = (sequence,chromosome,sequenceLength,Afrequency,Tfrequency,Gfrequency,Cfrequency, bimers, trimers, initiationrates)
arrayEsetA = []
arrayEsetT = []
arrayEsetG = []
arrayEsetC = []
# arrayEsetNothing =[]
sumEsetSequences = []

for j in range(maxLength):#max possible length of sequence, passing from every position in the sequence
    #print j
    As = Ts = Gs = Cs = 0
    countSkips = 0
    for i in range(topTen): #passing every sequence

        lenseq = int(sortedArray[i][2]) #length   !!!!PROSOXI!!!!!

        if j < lenseq: # if position j belongs to length
            if sortedArray[i][0][lenseq-(j+1)] == 'A':  #walk from right to left
                As += 1
            elif sortedArray[i][0][lenseq-(j+1)] == 'T':
                Ts += 1
            elif sortedArray[i][0][lenseq-(j+1)] == 'G':
                Gs += 1
            elif sortedArray[i][0][lenseq-(j+1)] == 'C':
                Cs += 1
        else:
            countSkips +=1 # Like this i count the sequences that dont have this position
        # print "sequence", sortedArray[i][0]
    pSetA = As/float(topTen)
    pSetT = Ts/float(topTen)
    pSetG = Gs/float(topTen)
    pSetC = Cs/float(topTen)
    pSetNothing = countSkips/float(topTen) #Now empty positions also have a Value
    # print "As", As
    # print "Ts", Ts
    # print "Gs", Gs
    # print "Cs", Cs
    # print "countSkips", countSkips
    # print "pSetA",pSetA
    # print "pSetT",pSetT
    # print "pSetG",pSetG
    # print "pSetC",pSetC
    # print "pSetNothing",pSetNothing


    Ab = Tb = Gb = Cb = Nothingb = 0
    count = 0
    for z in range(topTen,len(sortedArray)): #now check for the specific position in the rest 90%
        lenseq = int(sortedArray[z][2]) #length   !!!!PROSOXI!!!!!
        
        if j < lenseq:
            if sortedArray[z][0][lenseq-(j+1)] == 'A':  #walk from right to left
                Ab += 1
            elif sortedArray[z][0][lenseq-(j+1)] == 'T':
                Tb += 1
            elif sortedArray[z][0][lenseq-(j+1)] == 'G':
                Gb += 1
            elif sortedArray[z][0][lenseq-(j+1)] == 'C':
                Cb += 1
        else:
            Nothingb +=1

    pBackA = Ab/float(ninety)
    pBackT = Tb/float(ninety)
    pBackG = Gb/float(ninety)
    pBackC = Cb/float(ninety)
    pBackNothing = Nothingb/float(ninety)
    # print "Ab",Ab
    # print "Tb",Tb
    # print "Gb",Gb
    # print "Cb",Cb
    # print "Nothingb",Nothingb
    # print "pBackA",pBackA
    # print "pBackT",pBackT
    # print "pBackG",pBackG
    # print "pBackC",pBackC
    # print "pBackNothing",pBackNothing


    #Now for each position calculate the Eset of each base
    if pBackA != 0 and (pSetA/pBackA) != 0:
        esetA = pSetA*math.log(pSetA/pBackA)
    else:
        esetA = 0
    arrayEsetA.append(esetA)


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

    if pBackNothing != 0 and (pSetNothing/pBackNothing) !=0 :
        esetNothing = pSetNothing*math.log(pSetNothing/pBackNothing)
    else:
        esetNothing = 0
    arrayEsetNothing.append(esetNothing)

########################~~Summations~~########################################
#Now pass from every sequence and swap every base with its eset value
#Add them all together and get a number for every sequence. Use it as feature
#Take into account the spaces too.
for s in range(len(actualArray)):
    lenseq = int(actualArray[s][2]) # Length    !!!!PROSOXI!!!!!
    esetSequence = []
    for j in range(maxLength):

        if j < lenseq:
            if actualArray[s][0][j] == 'A':
                esetSequence.append(arrayEsetA[j])
            elif actualArray[s][0][j] == 'T':
                esetSequence.append(arrayEsetT[j])
            elif actualArray[s][0][j] == 'G':
                esetSequence.append(arrayEsetG[j])
            elif actualArray[s][0][j] == 'C':
                esetSequence.append(arrayEsetC[j])
        else:
            esetSequence.append(arrayEsetNothing[j])
    Summation = sum(esetSequence)
    # print "Summation "+`s`, Summation
    sumEsetSequences.append(Summation) #Append the Sum of the values of each Sequence-array.
    # print "esetSequence", esetSequence
# print "sumEsetSequences", len(sumEsetSequences)
sumEsetSequences = np.array(sumEsetSequences) # is an array that contains the Summation of the Eset values per sequence

np.save("Stress/sumEsetSequencesAlexey", sumEsetSequences)



'''
a = []
b = []
c = []

for i in range(noRows): # no or Rows in ActualArray = Existing Genes
    for j in range(len(stressedGenes)):
    
        if (actualArrayWithNames[i][0]==stressedGenes[j][0]) and (stressedGenes[j][2] != '-1'):
            
            a.append(actualArrayWithNames[i][2:-1])  # Get the first Features 
            b.append(stressedGenes[j][2])   # Get the Stress Class
            c.append(actualArrayWithNames[i][-1])    # Get the Initiation Rates

            #print "-->", existingGenes[i][0], stressedGenes[j][2], actualArrayWithNames[i][-1]
a = np.array(a).astype(float)
b = np.array(b).astype(float)
c = np.array(c).astype(float)
print "a", a.shape
print "b", b.shape
print "c", c.shape
'''



        # Pipeline([
        #     ('feature_selection', linear_model.Lasso(alpha=0.1, normalize=True)),
        #     ('regression', RandomForestRegressor(n_estimators=1000))])



#For 80/20 CV
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state=1)

MAE : mean_absolute_error(y_test, prediction)

GradientBoostingRegressor(n_estimators=500, max_depth=4, learning_rate=0.1, loss='huber', min_samples_leaf=3,random_state=0, verbose=1)

plot_partial_dependence(est, X_train, features, ....)