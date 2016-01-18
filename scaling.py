import numpy as np
from scipy.signal import argrelextrema
from matplotlib.pyplot import plot, show, errorbar
import matplotlib.pylab as plt
name = (raw_input("Enter the name of the 1st input file: "))
data1 = np.loadtxt(name)
data1 = data1.T
x = data1[0]
x = [i-500.0 for i in x]
y = data1[1]
plt.scatter(x,y)
plt.show()
minima = 0
i = 1
while(minima==0):
	check = y[501+i] - y[500+i]
	if (check>0):
		print i
		minima = 1
	else:
		i = i+1