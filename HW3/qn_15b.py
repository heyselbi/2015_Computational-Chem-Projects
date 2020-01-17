from scipy.constants import k, N_A
from math import *
from numpy import array,loadtxt, empty, linspace
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle

x=linspace(0,5)

import numpy 

px=numpy.exp(-x**2)

plot(x,px)    
xlabel("Distance between ligand and protein")
ylabel("Probability of the distance")
suptitle("Question 15 part b", fontsize=20)
show()
