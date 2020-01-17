from scipy.constants import k, N_A
from math import *
from numpy import array,loadtxt, empty, linspace
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle

T=linspace(1,3000, 3000)

import numpy 

p1=1
p2=numpy.exp(-3*4184/(T*N_A*k))
p3=3*numpy.exp(-3*4184/(T*N_A*k))

Q=p1+p2+p3

p1norm=p1/Q
p2norm=p2/Q
p3norm=p3/Q

dave=p1norm + p2norm*3 + p3norm*sqrt(5)

plot(T,dave)    
xlabel("Temperature in K")
ylabel("Average end-to-end distance")
suptitle("Question 11 part b", fontsize=20)
show()  
