from myfunctions import findBasesFrequency
a = 'CCTTCTTCTTTCTTAAAAAGTCTTAGTACGATTGACCAAGTCAGAAAAAAAAAAAAAAAGGAACTAAAAAAAGTTTTAATTAATTATGAGAGCTTTGGCATATTTCAAGAAGGGTGATATTCACTT'

Afrequency, Tfrequency, Gfrequency, Cfrequency = findBasesFrequency(a)
print Afrequency, Tfrequency, Gfrequency, Cfrequency
exit()
# print range(5)
# print ""
# for i in range(10):
# 	print "i", i
# 	if i in range(5):
# 		print "yes"
# 	print ""
# exit()
import numpy as np
sortedArray = np.load('Gritsenko/sortedArray.npy')
print "sortedArray  loaded, is the actualArrayWithNames Array (without the names) sorted by Init Rates: ", sortedArray.shape

temp = sortedArray[:,2].astype(int) # LEnGTH
maxLength=temp.max(axis=0)
print "maxLength", maxLength # Find max length
minLength=temp.min(axis=0)
print "minLength", minLength
print ""

arrayEsetA = []
arrayEsetT = []
arrayEsetG = []
arrayEsetC = []

#Pass from every position
for k in range(maxLength):
	print "position", k
	tempSequences = []

	#Select the sequences that have this length
	for s in range(len(sortedArray)):

		#Find Sequence's Length
		lenseq = int(sortedArray[s][2])
		print "lenseq", lenseq

		#If position belongs to the sequence select it
		if k in range(lenseq):
			tempSequences.append(sortedArray[s])
	tempSequences = np.array(tempSequences)
	print "tempSequences", tempSequences.shape
	print ""
	#Take the 10% of the sequences with the highest Initiation Rate
	topTen = int(len(tempSequences)*0.1)
	ninety = len(tempSequences) - topTen
	print "topTen",topTen
	print "ninety", ninety

	As = Ts = Gs = Cs = 0
	for i in range(topTen):
		lenseq = int(tempSequences[i][2])

		if tempSequences[i][0][lenseq-(k+1)] == 'A':  #walk from right to left
			As += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'T':
			Ts += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'G':
			Gs += 1
		elif tempSequences[i][0][lenseq-(k+1)] == 'C':
			Cs += 1
		else:
			print "error"

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
	exit()