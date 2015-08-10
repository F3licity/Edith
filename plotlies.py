# import plotly.plotly as py
# from plotly.graph_objs import *
# py.sign_in('Eftychia', '2puhmq6aj8')

# import matplotlib as m
# import matplotlib.pyplot as plt
# import numpy as np

# cdict = {
#   'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
#   'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
#   'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
# }

# cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

# x = np.arange(0, 10, .1)
# y = np.arange(0, 10, .1)
# X, Y = np.meshgrid(x,y)

# data = 2*( np.sin(X) + np.sin(3*Y) )

# def do_plot(n, f, title):
#     #plt.clf()
#     plt.subplot(1, 3, n)
#     plt.pcolor(X, Y, f(data), cmap=cm, vmin=-4, vmax=4)
#     plt.title(title)
#     plt.colorbar()

# plt.figure()
# do_plot(1, lambda x:x, "all")
# do_plot(2, lambda x:np.clip(x, -4, 0), "<0")
# do_plot(3, lambda x:np.clip(x, 0, 4), ">0")
# plt.show()

import matplotlib.pyplot as plt
# plt.figure(1)                # the first figure
# plt.subplot(211)             # the first subplot in the first figure
# plt.plot([1,2,3])
# plt.subplot(212)             # the second subplot in the first figure
# plt.plot([4,5,6])


# plt.figure(2)                # a second figure
# plt.plot([4,5,6])            # creates a subplot(111) by default

# plt.figure(1)                # figure 1 current; subplot(212) still current
# plt.subplot(211)             # make subplot(211) in figure1 current
# plt.title('Easy as 1,2,3')   # subplot 211 title

import matplotlib.pyplot as plt
# plot a line, implicitly creating a subplot(111)
plt.plot([1,2,3])

plt.subplot(211)
plt.plot(range(12))
plt.subplot(212, axisbg='y') # creates 2nd subplot with yellow background
plt.show()