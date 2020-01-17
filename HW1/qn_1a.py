from pylab import plot,show
from numpy import arange
from scipy.constants import k

#k is the Boltzmann constant, o is for Oxygen and cl is for Chlorine

sigmao=0.295e-9
sigmacl=0.335e-9
sigmaocl=(sigmao+sigmacl)/2
sigmaoo=(2*sigmao)/2
sigmaclcl=(2*sigmacl)/2

epso=61.6*k
epscl=173.5*k
epsocl=(epso*epscl)**0.5

#Potential Energy between Oxygen and Oxygen
Voo=[]
distanceoo=[]

for r in arange(1e-10,15e-10,1e-11):
    distanceoo.append(r)
    Voo.append(4*epso*((sigmaoo/r)**12-(sigmaoo/r)**6))
print(distanceoo)
print(Voo)
plot(distanceoo,Voo,"g-")

#Potential Energy between Chlorine and Chlorine
Vclcl=[]
distanceclcl=[]

for r in arange(1e-10,15e-10,1e-11):
    distanceclcl.append(r)
    Vclcl.append(4*epscl*((sigmaclcl/r)**12-(sigmaclcl/r)**6))
print(distanceclcl)
print(Vclcl)
plot(distanceclcl,Vclcl,"r--")


#Potential Energy between Oxygen and Chlorine
Vocl=[]
distanceocl=[]

for r in arange(1e-10,15e-10,1e-11):
    distanceocl.append(r)
    Vocl.append(4*epsocl*((sigmaocl/r)**12-(sigmaocl/r)**6))
print(distanceocl)
print(Vocl)
plot(distanceocl,Vocl,"ko")


show()
