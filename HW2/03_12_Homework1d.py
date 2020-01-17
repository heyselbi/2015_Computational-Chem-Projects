from math import *
from numpy import array,loadtxt, empty, linspace
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle
from scipy.constants import k, N_A

conf=loadtxt("N10.txt",int)

x=[0,4,10,16,18] #x coordinates of H proteins

totalE=[]

for j in range(0,4067):  #select row
    E=0
    for i in x: #select x coordinates of H proteins
        Hxo=conf[j,i]
        n=i+1 #y coordinates of H proteins
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

evenx=[0,2,4,6,8,10,12,14,16,18]

Rgall=[]

for a in range(0,4067):
    rijtotal=0.0
    for c in evenx: #x coordinate of aa
        xo=conf[a,c]
        f=c+1 #y coordinates of aa
        yo=conf[a,f]
        for d in evenx:
            if d-c>0:
                xf=conf[a,d]
                g=d+1
                yf=conf[a,g]
                rij=(xo-xf)**2 + (yo-yf)**2
                rijtotal+=rij
    Rgrow=rijtotal/100 #N=10, Rg of the row/conformation
    Rgall.append(Rgrow) #Rgall is the list of Rg's of all conformations
print "average gyration is: ", sum(Rgall)/len(Rgall)
    
for b in range(4067):
    if totalE[b]==0:
        Rgp0=(pnorm0*Rgall[b])/degeneracy[0]
    elif totalE[b]==-1:
        Rgp1=(pnorm1*Rgall[b])/degeneracy[1]
    elif totalE[b]==-2:
        Rgp2=(pnorm2*Rgall[b])/degeneracy[2]
    elif totalE[b]==-3:
        Rgp3=(pnorm3*Rgall[b])/degeneracy[3]
    else:
        Rgp4=(pnorm4*Rgall[b])/degeneracy[4]
        
Rgptotal=Rgp0+Rgp1+Rgp2+Rgp3+Rgp4
plot(T,Rgptotal)    
xlabel("Temperature in K")
ylabel("Radius of gyration")
suptitle("Radius of gyration change as a function of Temperature", fontsize=20)
show()  
  
