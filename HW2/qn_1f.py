from math import *
from numpy import array,loadtxt, empty, linspace
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle
from scipy.constants import k, N_A

conf=loadtxt("N10.txt",int)

x=[0,4,10,16,18] #x coordinates of H proteins

totalE=[]

for j in range(0,4067):  #select row
    E=0
    for d in x: #select x coordinates of H proteins
        Hxo=conf[j,d]
        n=d+1 #y coordinates of H proteins
        Hyo=conf[j,n]
        for l in x:
            Hxf=conf[j,l]
            m=l+1
            Hyf=conf[j,m]
            r=( (Hxo-Hxf)**2 + (Hyo-Hyf)**2 )**0.5
            if r==1:
                E+=(-1*r) #adding total energy
    totalE.append(E/2 +1) #need to div by 2 because r between 0,0 and 1,1 is the same as 1,1 and 0,0
    #need to add one, as there is a bond between last two H


zero,one,two,three,four=0,0,0,0,0

for s in range(4067):
    if totalE[s]==0:
        zero+=1
    elif totalE[s]==-1:
        one+=1
    elif totalE[s]==-2:
        two+=1
    elif totalE[s]==-3:
        three+=1
    else:
        four+=1
        
degeneracy=[zero, one, two, three, four]
energy=[0.0,-1.0,-2.0,-3.0,-4.0]

import numpy

T=linspace(100,500,401)
p0=degeneracy[0]*numpy.exp(-energy[0]*4184/(N_A*k*T))
p1=degeneracy[1]*numpy.exp(-energy[1]*4184/(N_A*k*T))
p2=degeneracy[2]*numpy.exp(-energy[2]*4184/(N_A*k*T))
p3=degeneracy[3]*numpy.exp(-energy[3]*4184/(N_A*k*T))
p4=degeneracy[4]*numpy.exp(-energy[4]*4184/(N_A*k*T))
Q=p0+p1+p2+p3+p4

pnorm0=p0/Q
pnorm1=p1/Q
pnorm2=p2/Q
pnorm3=p3/Q
pnorm4=p4/Q
#print pnorm0+pnorm1+pnorm2+pnorm3+pnorm4

Tf=T+1
p0Tf=degeneracy[0]*numpy.exp(-energy[0]*4184/(N_A*k*Tf))
p1Tf=degeneracy[1]*numpy.exp(-energy[1]*4184/(N_A*k*Tf))
p2Tf=degeneracy[2]*numpy.exp(-energy[2]*4184/(N_A*k*Tf))
p3Tf=degeneracy[3]*numpy.exp(-energy[3]*4184/(N_A*k*Tf))
p4Tf=degeneracy[4]*numpy.exp(-energy[4]*4184/(N_A*k*Tf))
QTf=p0Tf+p1Tf+p2Tf+p3Tf+p4Tf

dlnQdT=(numpy.log(QTf)-numpy.log(Q))/(Tf-T)

Utotala=k*T**2*dlnQdT*N_A


for b in range(4067):
    if totalE[b]==0:
        U0=pnorm0*totalE[b]*4184
    elif totalE[b]==-1:
        U1=pnorm1*totalE[b]*4184
    elif totalE[b]==-2:
        U2=pnorm2*totalE[b]*4184
    elif totalE[b]==-3:
        U3=pnorm3*totalE[b]*4184
    else:
        U4=pnorm4*totalE[b]*4184

Utotalb=U0+U1+U2+U3+U4

plot(T,Utotala,"r",label="Finite Difference Method")
plot(T,Utotalb,"b",label="Expectation Value Formula")    
xlabel("Temperature in K")
ylabel("Total Internal Energy per protein in Joules")
suptitle("Total Internal Energy as a function of Temperature", fontsize=20)
show()
