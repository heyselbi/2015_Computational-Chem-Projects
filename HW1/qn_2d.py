from math import *
from numpy import array
from pylab import plot,show,ylabel,xlabel

deltat=array([1e-4,1e-3,1e-2,1e-1,1],float)

pverlet=[]
s=0
xverlet=[]
xo=0
po=2
Uo=2
m=1
H=[0]
c=[0]
time=0
dHnew=0
dH=0

#deltat=1e-4
for k in range(10000000):
    a=float(xo+deltat[0]*(po/m)-(deltat[0]**2 * Uo * xo)/(2*m))
    xverlet.append(a)
    b=float(po-(deltat[0]*Uo*(xo+xverlet[s]))/2)
    pverlet.append(b)
    dH=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    xo=float(a)
    po=float(b)
    s+=1
    dHnew=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    deltaH=float((dHnew-dH)/deltat[0])
    H.append(deltaH)
    time+=deltat[0]
    c.append(time)
plot(c,H)

#deltat=1e-3
pverlet=[]
s=0
xverlet=[]
xo=0
po=2
Uo=2
m=1
H=[0]
c=[0]
time=0
dHnew=0
dH=0
for k in range(1000000):
    a=float(xo+deltat[1]*(po/m)-(deltat[1]**2 * Uo * xo)/(2*m))
    xverlet.append(a)
    b=float(po-(deltat[1]*Uo*(xo+xverlet[s]))/2)
    pverlet.append(b)
    dH=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    xo=float(a)
    po=float(b)
    s+=1
    dHnew=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    deltaH=float((dHnew-dH)/deltat[1])
    H.append(deltaH)
    time+=deltat[1]
    c.append(time)
plot(c,H)

#deltat=1e-2
pverlet=[]
s=0
xverlet=[]
xo=0
po=2
Uo=2
m=1
H=[0]
c=[0]
time=0
dHnew=0
dH=0
for k in range(100000):
    a=float(xo+deltat[2]*(po/m)-(deltat[2]**2 * Uo * xo)/(2*m))
    xverlet.append(a)
    b=float(po-(deltat[2]*Uo*(xo+xverlet[s]))/2)
    pverlet.append(b)
    dH=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    xo=float(a)
    po=float(b)
    s+=1
    dHnew=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    deltaH=float((dHnew-dH)/deltat[2])
    H.append(deltaH)
    time+=deltat[2]
    c.append(time)
plot(c,H)

#deltat=1e-1
pverlet=[]
s=0
xverlet=[]
xo=0
po=2
Uo=2
m=1
H=[0]
c=[0]
time=0
dHnew=0
dH=0
for k in range(10000):
    a=float(xo+deltat[3]*(po/m)-(deltat[3]**2 * Uo * xo)/(2*m))
    xverlet.append(a)
    b=float(po-(deltat[3]*Uo*(xo+xverlet[s]))/2)
    pverlet.append(b)
    dH=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    xo=float(a)
    po=float(b)
    s+=1
    dHnew=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    deltaH=float((dHnew-dH)/deltat[3])
    H.append(deltaH)
    time+=deltat[3]
    c.append(time)
plot(c,H)

#deltat=1
pverlet=[]
s=0
xverlet=[]
xo=0
po=2
Uo=2
m=1
H=[0]
c=[0]
time=0
dHnew=float(0)
dH=float(0)
for k in range(1000):
    a=float(xo+deltat[4]*(po/m)-(deltat[4]**2 * Uo * xo)/(2*m))
    xverlet.append(a)
    b=float(po-(deltat[4]*Uo*(xo+xverlet[s]))/2)
    pverlet.append(b)
    dH=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    xo=float(a)
    po=float(b)
    s+=1
    dHnew=float(po**2/(2*m) + (Uo*(xo**2 -2))/2)
    deltaH=float((dHnew-dH)/deltat[4])
    H.append(deltaH)
    time+=deltat[4]
    c.append(time)
plot(c,H)

xlabel("Time")
ylabel("dH/dt")
show()

#the optimum step size to use is the one that is the smallest (ie. 1e-4)
#it is because the smaller increments of change in Hamiltonian are recorded, while bigger step size are less accurate
#it is obvious from the plot, the 1e-4 curve "oscillates" in smaller range of Hamiltonian values --> more accurate and optimum
