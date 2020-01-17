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
    totalE.append(E/2 +1) #need to div 2 because r between 0,0 and 1,1 is the same as 1,1 and 0,0
    #need to add one, as there is a bond between last two H


from collections import defaultdict
d = defaultdict(list)
for i, x in enumerate(totalE):
    d[x].append(i)
    
k = min(d.keys())
minErow=array(d[k])

print "Minimum energy observed is ", min(totalE), "eps kcal/mol"
print "Row number of minimum E observed: ", minErow+1 #Because 0 is also number, 1 is added

a=minErow.size

for q in range(a):
    minimum=minErow[q]
    xrow=[]
    yrow=[]
    for n in range(0,20):
        if n%2==0:
            xrow.append(conf[minimum,n])
            p=n+1
            yrow.append(conf[minimum,p])

    #print(xrow,yrow)
    plot(xrow,yrow)
ylim(-4,4)
xlim(-1,3)
xlabel("x coordinate", fontsize=16)
ylabel("y coordinate", fontsize=16)
suptitle("Protein shapes of the ones showing lowest energy state", fontsize=20)
show()
