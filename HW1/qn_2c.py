from math import *
from numpy import array

deltat=[1e-4,1e-3,1e-2,1e-1,1]

xverlet=[] #list of positions at deltat
xo=0
po=2
Uo=2
m=1

for i in deltat:
    a=xo+i*(po/m)+(i**2 * Uo * xo)/(2*m)
    xverlet.append(a)
print(xverlet)


pverlet=[]
s=0

for j in deltat:
    b=po-(j*Uo*(xo+xverlet[s]))/2
    s+=1
    pverlet.append(b)
print(pverlet)
