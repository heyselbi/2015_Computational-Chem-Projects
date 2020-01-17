from math import *
from numpy import array,loadtxt
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle

conf=loadtxt("N10.txt",int)

x=[0,4,10,16,18] #x coordinates of H proteins

totalE=[]

for j in range(0,4067):  #select row
    E=0
    for i in x: #select x coordinates of H proteins
        Hxo=conf[j,i]
        k=i+1 #y coordinates of H proteins
        Hyo=conf[j,k]
        for l in x:
            Hxf=conf[j,l]
            m=l+1
            Hyf=conf[j,m]
            r=( (Hxo-Hxf)**2 + (Hyo-Hyf)**2 )**0.5
            if r==1:
                E+=(-1*r) #adding total energy
    totalE.append(E/2 +1) #need to div by 2 because r between 0,0 and 1,1 is the same as 1,1 and 0,0
    #need to add one, as there is a bond between last two H


b=array(totalE)

zero,one,two,three,four=0,0,0,0,0

for s in range(4067):
    if b[s]==0:
        zero+=1
    elif b[s]==-1:
        one+=1
    elif b[s]==-2:
        two+=1
    elif b[s]==-3:
        three+=1
    else:
        four+=1
        
degeneracy=array([zero, one, two, three, four],int)
energy=array([0,-1,-2,-3,-4],int)

plot(energy,degeneracy)
ylim(-1,3000)
xlim(-4.5,0.5)
xlabel("Energy in epsilon", fontsize=16)
ylabel("Degeneracy", fontsize=16)
suptitle("Energy Vs Degeneracy of Protein Conformations", fontsize=20)
show()
