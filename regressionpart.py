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
    allcChromosomes = featureArray[featureArray[:,0] == float(c)] # np.where(featureArray[:,1]==0) Take all the rows per chromosome

    X = allcChromosomes[:,1:-1] #Features
    y = allcChromosomes[:,-1] #initiation rates
    # print "X", X.shape
    # print "y", y.shape


    ##### Calculate the Correlation among the features and the target

    pearsonsCorrelations = []
    spearmanCorrelations = []
    # transposedX = X.T
    if len(y) != 0:
        for i in range(0, len(X.T)):
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
        # plt.show()
        plt.savefig("CorrelationCoefficientFor" + `c`+"Chrom"+ `csvFile` + ".png")
        plt.clf()


        # Do SVR
        skf = cross_validation.KFold(len(y),n_folds=5)
        allScores = []
        print "X and Y", X.shape, y.shape
        fold = 1
        for train_index, test_index in skf:
            #print("TRAIN:", train_index, "TEST:", test_index)
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            svr_linear = svm.SVR(kernel='linear', degree=4, C=1, cache_size=300) #SupportVectorMachine - Support Vector Regression  #RandomForestRegressor(n_estimators=100) 
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


            #Images
            plt.close('all')
            plt.figure()
            plt.subplots_adjust(hspace=.5, wspace=.5)
            plt.subplot(1,2,1)
            plt.title('SVR train')
            allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=1, label='train')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_linear_predictTrain, c='r', edgecolor='none', s=1, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')

            plt.subplot(1,2,2)
            allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=1, label='test')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_linear_predictTest, c='r', edgecolor='none', s=1, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')
            plt.title('SVR test')
            plt.legend()
            plt.suptitle('SVR_AllFeats'+ `csvFile` +`fold`)
            plt.savefig('correctciandrinidata/SVR_AllFeatsFor' + `c` + 'Chrom' + `csvFile` +`fold`+'.svg')
            # plt.show()
            plt.clf()
            fold += 1


        # Do Random Forest
        allScores = np.array(allScores)
        print "SVR Scores" + `c`, allScores
        print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
        print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(svr_linear.predict(X_train), y_train)))
        print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(svr_linear.predict(X_test), y_test)))


        allScores = []
        fold = 1
        for train_index, test_index in skf:
            #print("TRAIN:", train_index, "TEST:", test_index)
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            rndmForest = RandomForestRegressor(n_estimators=100)
            print "Start training Random Forest"
            start_time_training = time.time()
            rndmForest.fit(X_train, y_train)
            print "Done training in ",time.time() - start_time_training, "seconds"
            score = rndmForest.score(X_test, y_test)
            allScores.append(score)
            y_predictTest = rndmForest.predict(X_test)
            y_predictTrain = rndmForest.predict(X_train)

            # #Random Forest
            # res = RandomForestRegressor(n_estimators=100)
            # res = res.fit(X_train, y_train)
            Y_predicted = rndmForest.predict(X_train)
            Y_test_predicted = rndmForest.predict(X_test)



           # print "y_test", y_test
            #print "y_predictTrain", y_predictTrain
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
            #exit()
            #Images
            plt.close('all')
            plt.figure()
            plt.subplots_adjust(hspace=.5, wspace=.5)
            plt.subplot(1,2,1)
            plt.title('rndmForest train scaled')
            allSampleIndexes=np.arange(len(X_train[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_train, c='g', edgecolor='none', s=1, label='train')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_predictTrain, c='r', edgecolor='none', s=1, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')

            plt.subplot(1,2,2)
            allSampleIndexes=np.arange(len(X_test[:])) #we put the len so to give to each sample a number
            plt.scatter(allSampleIndexes, y_test, c='g', edgecolor='none', s=1, label='test')
            plt.hold('on')
            plt.scatter(allSampleIndexes, y_predictTest, c='r', edgecolor='none', s=1, label='prediction')
            plt.xlabel('#of rows')
            plt.ylabel('initiation rates (log)')
            plt.title('rndmForest test')
            plt.legend()
            plt.suptitle('rndmForest_AllFeats'+ `csvFile` +`fold`)
            plt.savefig('correctciandrinidata/rndmForest_AllFeatsFor' + `c` + 'Chrom' + `csvFile` +`fold`+'.svg')
            # plt.show()
            plt.clf()
            fold += 1

        # print "Prediction of X is ", Y_predicted
        allScores = np.array(allScores)
        print "Random Forest Scores" + `c`, allScores
        print("Accuracy: %0.2f (+/- %0.2f)" % (allScores.mean(), allScores.std()))
        print('MeanAbsoluteError Train: {}'.format(metrics.mean_absolute_error(rndmForest.predict(X_train), y_train)))
        print('MeanAbsoluteError Test: {}'.format(metrics.mean_absolute_error(rndmForest.predict(X_test), y_test)))