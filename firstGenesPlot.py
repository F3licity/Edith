import numpy as np
import matplotlib.pyplot as plt
import csv

def graph():
	genes = []
	InitiationRates = []
	with open('eftychia.csv', 'rb') as csvfile:
		genereader = csv.reader(csvfile, delimiter=',')

		for row in genereader:
			print ', '.join(row)
			genes.append(row[0])
			InitiationRates.append(float(row[1]))

		# num_above = 0.0001
		# for row in genereader:
		# 	num_above = [val+1 if key <= i else val for key,val in enumerate(num_above)]
		# 	print ', '.join(num_above)
		# 	genes.append(num_above[0])
		# 	InitiationRates.append(float(num_above[1]))


	print "MEAN:",np.mean(InitiationRates)
	print "STD:",np.std(InitiationRates)
	plt.plot(InitiationRates,'ro')
	plt.show()

graph()


# for i in B:
#     num_above = [val+1 if key <= i else val for key,val in enumerate(num_above)]
