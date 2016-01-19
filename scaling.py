import numpy as np
from scipy.signal import argrelextrema
from matplotlib.pyplot import plot, show, errorbar
import matplotlib.pylab as plt
times = [10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0,4000.0,6000.0]
# times = [50.0,100.0,200.0,500.0,1000.0,2000.0,4000.0,6000.0]
# times = [200.0,500.0,1000.0,2000.0,4000.0,6000.0]
# times = [1000.0]
t = len(times)
x_minimas = np.zeros([t])
x_maximas = np.zeros([t])
peak_heights_heat = np.zeros([t])
peak_heights_trav = np.zeros([t])
f = open("maximas.dat", "w")
for i in xrange(t):
	name = "den-and-stddev-r0.10time"+ str(int(times[i])) + "avg-over1000.dat"
	data1 = np.loadtxt(name)
	# f = open('den-and-stddev-r0.10time1000avg-over1000.dat', 'r')
	# data1 = f.read()
	data1 = data1.T
	x = data1[0]
	x = [j-500.0 for j in x]
	y = data1[1]
	peak_heights_heat[i] = y[501]
	# print "peak_heights_heat: ", y[501]
	x_minima = 0
	x_maxima = 0
	j = 1
	while (x_minima==0 or x_maxima==0):
		check = y[501+j] - y[500+j]
		if (check>0 and x_minima==0):
			x_minimas[i] = j
			# print "x_minima: ", j
			x_minima = 1
			j = j+1
		elif (check<0 and x_minima==1):
			x_maximas[i] = j
			# print "x_maxima: ", j
			peak_heights_trav[i] = y[500+j]
			# print "peak_heights_trav: ", y[500+j]
			x_maxima = 1
		else:
			j = j+1
	f.write(str(times[i]) + '\t' + str(x_maximas[i]) + '\n')
f.close()		
	# plt.scatter(x,y)
	# plt.show()
# plt.scatter(times,x_minimas)
# plt.show()
# plt.scatter(times,x_maximas)
# plt.show()
print "\n" "Times considered: ", times, "\n"
slope_xmin, intercept_xmin = np.polyfit(np.log(times), np.log(x_minimas), 1)
print "\n", "X position of the minima scales with time as t to the power" , slope_xmin
slope_xmax, intercept_xmax = np.polyfit(np.log(times), np.log(x_maximas), 1)
print "X position of the travelling mode scales with time as t to the power" , slope_xmax
slope_heat, intercept_heat = np.polyfit(np.log(times), np.log(peak_heights_heat), 1)
print "Height of the heat mode peak scales with time as t to the power" , slope_heat
slope_trav, intercept_trav = np.polyfit(np.log(times), np.log(peak_heights_trav), 1)
print "Height of the travelling mode peak scales with time as t to the power" , slope_trav, "\n"