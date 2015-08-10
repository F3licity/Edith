import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import scipy
from scipy.stats.stats import pearsonr
from myfunctions import readCSV

csvFile = 'Ciandrini/YeastMineResult.csv'
delimiter = ','
goTermsCiandrini = readCSV(csvFile, delimiter)
goTermsCiandrini = np.array(goTermsCiandrini)
print goTermsCiandrini.shape
print goTermsCiandrini[0:5]

existingGenes = np.load('Ciandrini/existingGenes.npy')
# print "existingGenes", existingGenes.shape, existingGenes[0]
print ""
actualArray   = np.load('Ciandrini/actualArray.npy')
# print "actualArray",   actualArray.shape, actualArray[0]


#Make a copy array of the existingGenes(names) + a column of zeros
a = []
for i in range(len(existingGenes)):
	a.append((existingGenes[i][0], 0))
a = np.array(a)

#Check the name and replace the 0 with a value
#value = number of occurrences of the Name in the Stress Array
#in other words, with how many GO Terms (rel. STress)

for i in range(len(a)):
