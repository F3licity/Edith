# from pylab import *

# t = arange(0.0, 2.0, 0.01)
# s = sin(2*pi*t)
# plot(t, s)   #aksones sto kartesiano

# xlabel('time (s)')
# ylabel('voltage (mV)')
# title('About as simple as it gets, folks')
# grid(True)
# savefig("test.png")
# show()
####################################################
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# def update_line(num, data, line):
#     line.set_data(data[...,:num])
#     return line,

# fig1 = plt.figure()

# data = np.random.rand(2, 25)
# l, = plt.plot([], [], 'r-')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.xlabel('x')
# plt.title('test')
# line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#     interval=50, blit=True)
# #line_ani.save('lines.mp4')

# fig2 = plt.figure()

# x = np.arange(-9, 10)
# y = np.arange(-9, 10).reshape(-1, 1)
# base = np.hypot(x, y)
# ims = []
# for add in np.arange(15):
#     ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))

# im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
#     blit=True)
# #im_ani.save('im.mp4', metadata={'artist':'Guido'})

# plt.show()
##################################################################
# from StringIO import StringIO
# import numpy as np

# s = StringIO("1,1.3,abcde")
# s.seek(0) # needed for StringIO example only
# data = np.genfromtxt(s, dtype=None,
# 	names = ['myint','myfloat','mystring'], delimiter=",")
# print data, np.result_type(data)
##################################################
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


	#genes = np.genfromtxt('eftychia.csv', delimiter=',', names=['Genes', 'InitiationRates'])
	#genes = np.unique(genes)
	#genes, value = np.loadtxt('eftychia.csv', delimiter=',', unpack=True)

	print "MEAN:",np.mean(InitiationRates)
	print "STD:",np.std(InitiationRates)
	plt.plot(InitiationRates,'ro')
	plt.show()
	#fig = plt.figure()

	#from mpl_toolkits.mplot3d import axes3d
	#ax1 = fig.add_subplot(1,1,1, projection='3d')
	#plt.plot(genes, InitiationRates, '-')
	#plt.title('title')
	#Xuniques, X = np.unique(df['X'], return_inverse=True)
	#plt.ylabel('value')
	#plt.xlabel('genes')
	#plt.show()

graph()
