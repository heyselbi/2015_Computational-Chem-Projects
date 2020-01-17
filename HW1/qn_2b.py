from math import *
from numpy import array,empty

deltat=[1e-4,1e-3,1e-2,1e-1,1]


#SIMULATED
xsim=[] #list of positions at deltat
for i in deltat:
    a=2*i
    xsim.append(a)
print(xsim)


#EXACT
xexct=[]
for j in deltat:
    b=2*2**-0.5 * sin(j*2**0.5)
    xexct.append(b)
print(xexct)


#ERROR
error=[]
s=0
for k in deltat:
    c=(xsim[s]-xexct[s])/k**2
    error.append(c)
    s+=1
print(error)


#PLOTTING
from pylab import plot, show, xlim, ylim
plot(error,deltat,"g")
xlim(-0.5,1)
ylim(-0.5,1.5)
show()
