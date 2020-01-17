from numpy import loadtxt,array
from math import *

values=loadtxt("O2CL2.txt",float)  #2nd column 1 is O2 and 2 is CL

number=values[:,0]
atom=values[:,1]
x=values[:,2]
y=values[:,3]
z=values[:,4]

from scipy.constants import k

sigmao=0.295
sigmacl=0.335
sigmaocl=(sigmao+sigmacl)/2

epso=61.6*k
epscl=173.5*k
epsocl=(epso*epscl)**0.5

#7 is Chlorine, 1 is Ox and 2 Cl
Fxseven,Fyseven,Fzseven=0,0,0

for number in range(200):
    if atom[number]==1 and number!=6:
        r=( (x[number]-values[6,2])**2 + (y[number]-values[6,3])**2 + (z[number]-values[6,4])**2 )**0.5
        Fx=24*epsocl*(x[number]-values[6,2])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fxseven+=Fx
        Fy=24*epsocl*(y[number]-values[6,3])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fyseven+=Fy
        Fz=24*epsocl*(z[number]-values[6,4])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fzseven+=Fz
        
    elif atom[number]==2 and number!=6:
        r=((x[number]-values[6,2])**2 + (y[number]-values[6,3])**2 + (z[number]-values[6,4])**2 )**0.5
        Fx=24*epscl*(x[number]-values[6,2])*((20*(sigmacl**12)/r**14) - (10*(sigmacl**6)/r**8))
        Fxseven+=Fx
        Fy=24*epscl*(y[number]-values[6,3])*((20*(sigmacl**12)/r**14) - (10*(sigmacl**6)/r**8))
        Fyseven+=Fy
        Fz=24*epscl*(z[number]-values[6,4])*((20*(sigmacl**12)/r**14) - (10*(sigmacl**6)/r**8))
        Fzseven+=Fz
        
    else:
        Fxseven+=0
        Fyseven+=0
        Fzseven+=0
        
Fseven=[Fxseven, Fyseven, Fzseven]

Fsevenmag= (Fxseven**2 + Fyseven**2 + Fzseven**2 )**0.5

print("The total force on atom 7 is: ", Fseven, "Newtons", "and magnitude is: ", Fsevenmag, "Newtons")


#25 is Oxygen, 1 is Ox and 2 Cl

Fxtf,Fytf,Fztf=0,0,0

for number in range(200):
    if atom[number]==1 and number!=24:
        r=( (x[number]-values[24,2])**2 + (y[number]-values[24,3])**2 + (z[number]-values[24,4])**2 )**0.5
        Fx=24*epso*(x[number]-values[24,2])*((20*(sigmao**12)/r**14) - (10*(sigmao**6)/r**8))
        Fxtf+=Fx
        Fy=24*epso*(y[number]-values[24,3])*((20*(sigmao**12)/r**14) - (10*(sigmao**6)/r**8))
        Fytf+=Fy
        Fz=24*epso*(z[number]-values[24,4])*((20*(sigmao**12)/r**14) - (10*(sigmao**6)/r**8))
        Fztf+=Fz
        
    elif atom[number]==2 and number!=24:
        r=((x[number]-values[24,2])**2 + (y[number]-values[24,3])**2 + (z[number]-values[24,4])**2 )**0.5
        Fx=24*epsocl*(x[number]-values[24,2])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fxtf+=Fx
        Fy=24*epsocl*(y[number]-values[24,3])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fytf+=Fy
        Fz=24*epsocl*(z[number]-values[24,4])*((20*(sigmaocl**12)/r**14) - (10*(sigmaocl**6)/r**8))
        Fztf+=Fz
        
    else:
        Fxtf+=0
        Fytf+=0
        Fztf+=0
        
Ftf=[Fxtf, Fytf, Fztf]
Ftfmag= (Fxtf**2 + Fytf**2 + Fztf**2 )**0.5
print("The total force on atom 25 is: ", Ftf, "Newtons", "and magnitude is: ", Ftfmag, "Newtons")
