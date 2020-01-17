from math import *
from numpy import array,loadtxt, empty, linspace
from pylab import plot,show,ylabel,xlabel, ylim, xlim, suptitle
from scipy.constants import k, N_A
import itertools

conf=loadtxt("N10.txt",int)

randomx=list(itertools.combinations([0,2,4,6,8,10,12,14,16,18],5))

talllist=[]
Cvalllist=[]
for b in range(len(randomx)):
    x=randomx[b] #x coordinates of H proteins

    totalE=[]

    for j in range(0,4067):  #select row
        E=0
        for d in x: #select x coordinates of H proteins
            Hxo=conf[j,d]
            n=d+1 #y coordinates of H proteins
            Hyo=conf[j,n]
            for l in x:
                diff=0
                diff=fabs(d-l)
                if diff>2:
                    Hxf=conf[j,l]
                    m=l+1
                    Hyf=conf[j,m]
                    r=( (Hxo-Hxf)**2 + (Hyo-Hyf)**2 )**0.5
                    if r==1:
                        E+=(-1*r) #adding total energy
        totalE.append(E/2) #need to div by 2 because r between 0,0 and 1,1 is the same as 1,1 and 0,0
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
    
    Utotal=k*Tf**2*dlnQdT*N_A
    
    Cvall=[]
    h=1
    for g in range(400):
        Cv=Utotal[h]-Utotal[g]
        Cvall.append(Cv)
        h+=1
    
    #plot(T[1:401],Cvall,"r")
    
    #print "Maximum value of heat capacity/Melting point is: ", max(Cvall), "Joules"
    Cvalllist.append(max(Cvall))
    maxCvall=max(Cvall)
    indexmaxCvall= [i for i, j in enumerate(Cvall) if j == maxCvall]
    t=T[indexmaxCvall[0]]
    talllist.append(t)
            
#xlabel("Temperature in Kelvin")
#ylabel("Heat capacity in Joules")
#suptitle("Heat capacity as a function of Temperature", fontsize=20)
#show()

print "highest melting T detected: ", max(talllist)
print "highest Cv detected: ", max(Cvalllist)

maxCv = max(Cvalllist)
indexmax= [i for i, j in enumerate(Cvalllist) if j == maxCv]
#print indexmax

lenindexmax=len(indexmax)
for a in range(lenindexmax): 
    p=indexmax[a]
    print "x coordinates of possible configurations of H aa's: ", randomx[p]
    
maxt = max(talllist)
indexmaxt= [i for i, j in enumerate(talllist) if j == maxt]


lenindexmaxt=len(indexmaxt)
print lenindexmaxt
for c in range(lenindexmaxt): 
    q=indexmaxt[c]
    print "x coordinates of possible configurations of H aa's: ", randomx[q]
