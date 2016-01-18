import numpy as np
from scipy.signal import argrelextrema
from matplotlib.pyplot import plot, show, errorbar
import matplotlib.pylab as plt
times = [10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0,4000.0,6000.0]
t = len(times)
minimas = np.zeros([t])
for i in xrange(t):
	name = "den-and-stddev-r0.10time"+ str(int(times[i])) + "avg-over1000.dat"
	data1 = np.loadtxt(name)
	# f = open('den-and-stddev-r0.10time1000avg-over1000.dat', 'r')
	# data1 = f.read()
	data1 = data1.T
	x = data1[0]
	x = [j-500.0 for j in x]
	y = data1[1]
	# plt.scatter(x,y)
	# plt.show()
	minima = 0
	j = 1
	while(minima==0):
		check = y[501+j] - y[500+j]
		if (check>0):
			minimas[i] = j
			minima = 1
		else:
			j = j+1
plt.scatter(times,minimas)
plt.show()

slope, intercept = np.polyfit(np.log(times), np.log(minimas), 1)
print(slope)