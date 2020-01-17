from math import *
from numpy import array,loadtxt, empty, linspace, exp
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle

k=0.0000031668114 #Boltzmann in Hartree/K
angle=linspace(0,360,37) #dihedral angles from 0 to 360
energy=loadtxt("04_24_Homework43c.txt",float) #all values are in Hartree

#########################
#Homework 4 - Question 3c
"""plot(angle, energy, "r-")
xlabel("Dihedral angle in Degrees")
ylabel("Energy in Hartree")
suptitle("Energy versus Dihedral angle", fontsize=20)
xlim(-5,365)
show()"""

#########################
#Question 3d
"""T=(270,300,400)
import numpy
for t in range(3): #choose temperature
    Q=0.0
    prob=[]
    for i in range(37): #calculate Q
        pq=0
        pq=numpy.exp(-energy[i]/(1000*k*T[t]))
        Q+=pq
    for j in range(37): #calculate normalized probabilities
        p=0
        p=numpy.exp(-energy[j]/(1000*k*T[t]))/Q
        prob.append(p)
    plot(angle,prob)
    plot.label=T[t], "K"
    
xlim(-5,365)   
xlabel("Dihedral angle in Degrees")
ylabel("Normalized probability")
suptitle("Populations versus dihedral angle", fontsize=20)
show()"""

###########################
#Question 3e
T=linspace(270,400,131)
ptotal=[]
import numpy
for i in range(131): 
    Q=0.0   
    thetatotal=0.0
    for j in range(37): #calculate Q
        pq=0
        pq=numpy.exp(-energy[j]/(1000*k*T[i]))
        Q+=pq
    for a in range(37): #calculate normalized probabilities
        p=0
        p=angle[a]*numpy.exp(-energy[a]/(1000*k*T[i]))/Q
        thetatotal+=p
    ptotal.append(thetatotal)
plot(T,ptotal)
xlabel("Temperature in K")
ylabel("Average torsional angle in degrees")
suptitle("Torsional angle versus Temperature", fontsize=20)
show()
