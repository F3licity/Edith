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
print "goTermsCiandrini",goTermsCiandrini.shape
#['GO:0000001' 'YAL048C']


existingGenes = np.load('Ciandrini/existingGenes.npy')
# print "existingGenes", existingGenes.shape, existingGenes[0]
print ""
actualArray   = np.load('Ciandrini/actualArray.npy')
# print "actualArray",   actualArray.shape, actualArray[0]


#Make a copy array of the existingGenes(names) + a column of zeros
genes = []
for i in range(len(existingGenes)):
	genes.append((existingGenes[i][0], 0))
genes = np.array(genes)
#['YAL062W' '0']


terms = np.load('stressALLGoTerms.npy')
#['GO:0000012' 'GO:0000045' 'GO:0000046' 'GO:0000077' 'GO:0000161' etc]


#Check the name and replace the 0 with a value
#value = number of occurrences of the Name in the Stress Array
#in other words, how many GO Terms (rel. STress) a gene has.
#Loop through the array Of Genes and Goterms as they come from YeastMineQuery
for i in range(len(goTermsCiandrini)):

	#If the GoTerm the Gene has is related to Stress
	if goTermsCiandrini[i][0] in terms:

		#Find the Position of the gene in the array and increase the value
		index = np.where(genes[:,0] == goTermsCiandrini[i][1])
		index = index[0][0]
		# print index
		# print genes[index][1]
		genes[index][1] = genes[index][1].astype(int)+1

print genes
np.save('Ciandrini/stressedGenes', genes)

temp = np.array(genes[:,1],dtype=int)
maxValue = temp.max(axis=0) #axis=0:column: Indicates the no of Points/rows.
minValue = temp.min(axis=0)
print "Assigned Stress Terms Per gene"
print "maxValue", maxValue
print "minValue", minValue

