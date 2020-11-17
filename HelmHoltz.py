import numpy as np 
import matplotlib.pyplot as plt
from scipy import special as sp 

def singlecoil(a, r):
    "Magnetic field due to a single coil"
    # The coil is shifted from origin by a distance d along z-axis
    # Hence, the center of the coil is at (0,0,d) with normal along +ve z-axis
    # a is the radius of the coil in units of d
    # r is the list of coordinates where magnetic field is to be calculated in units of d
    # r is in general of form (x,y,z) w.r.t origin in units of d

    x=r[0]
    y=r[1]
    z=r[2]
    rho=np.sqrt(x**2+y**2)
    B=[0, 0]
    k=np.sqrt((4*a*rho)/((a+rho)**2+(z-1)**2))
    K=sp.ellipk(k**2)
    E=sp.ellipe(k**2)
    B[0]=(1/np.sqrt((a+rho)**2+(z-1)**2))*(K+E*(a**2-rho**2-(z-1)**2)/((a-rho)**2+(z-1)**2))
    if rho==0:
        B[1]=0
    else:
        B[1]=((z-1)/rho)*(1/np.sqrt((a+rho)**2+(z-1)**2))*(-K+E*(a**2-rho**2-(z-1)**2)/((a-rho)**2+(z-1)**2))
    return B

def antihelm(a, r):
    "Anti-Helmholtz Coils"
    Bres=[0, 0]
    x=r[0]
    y=r[1]
    z=r[2]
    Bres1=singlecoil(a, [x, y, z])
    Bres2=singlecoil(a, [x, y, -z])
    Bres[0]=Bres1[0]-Bres2[0]
    Bres[1]=Bres1[1]-Bres2[1]
    return Bres 

def helm(a, r):
    "Helmholtz Coils"
    Bres=[0, 0]
    x=r[0]
    y=r[1]
    z=r[2]
    Bres1=singlecoil(a, [x, y, z])
    Bres2=singlecoil(a, [x, y, -z])
    Bres[0]=Bres1[0]+Bres2[0]
    Bres[1]=Bres1[1]+Bres2[1]
    return Bres 

def appgrad(a, d):
    "Approximate Constant Gradient"
    Cgrad=6*np.pi*a/(a**2+1)**2.5/d
    return Cgrad

def linearField(a, r):
    "Approximate Linear Field"
    linB=6*np.pi*a*r[2]/(a**2+1)**2.5
    return linB


#Call the function

def adjustB(a, z, d):
    "Ajusted Field"
    #z = np.linspace(-1,1,1000)
    Badj=2*10**(-5)*antihelm(a, [0, 0, z])[0]/d
    return Badj

def adjustB2(a, z, d):
    "Ajusted Field"
    #z = np.linspace(-1,1,1000)
    Badj=2*10**(-3)*helm(a, [0, 0, z])[0]/d
    return Badj

Z=np.linspace(-0.04,0.04,50000)
Bfin=np.zeros(50000)
len=0.0
for i in range(0,24):
    d=0.002*i+0.084
    for j in range(0,17):
        rad=0.037+0.002*j
        a=rad/d
        z2=Z/d
        len=len+2*np.pi*rad
        Bfin=Bfin+adjustB2(a, z2, d)
print(len)
#for i in range(0,20):
#    d=0.002*i+0.035
#    for j in range(0,20):
#        rad=0.035+0.002*j
#        a=rad/d
#        z2=Z/d
#        len=len+2*np.pi*rad
#        Bfin=Bfin+adjustB2(a, z2, d)
#print(len)
z = np.linspace(-1,1,1000)
Bnew=adjustB(0.71, z,  0.045)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_xlabel('Distance from center of system (meters)')
ax.set_ylabel('Magnetic Field in Gauss/A')
ax.set_ylim(0.0, 20.0)

#for i in range(0,1500):
#    a0=i/500+0.25
#    Bplot=adjustB(a0, z, 0.05)
#    plt.plot(z, np.gradient(Bplot,z), 'r')
# plot the function
#plt.plot(Z, Bfin, 'r')
#plt.plot(z, Bnew, 'r')
plt.plot(Z, Bfin, 'k')
print (Bfin[25001])

"""
B=antihelm(0.5, [0, 0, 0])
print (B)   

z = np.linspace(-1,1,1000)
R = np.linspace(0.25,2.5,1000)
#Bz=antihelm(R,[0, 0, 1])
#print R, Bz
Bz=np.zeros([1000, 1000])
grad=np.zeros(1000)
max=-np.inf
dgrad=np.zeros([1000, 1000])
k=0
for i in range(0,999):
    r0=i/500+0.25
    Bz[i]=antihelm(r0,[0, 0, z])[0]
    grad[i]=np.gradient(Bz[i])[500]
    dgrad[i]=np.gradient(np.gradient(Bz[i]))
    if grad[i]>max:
        max=grad[i]
        k=i
print (max, k/100+0.25, grad[100], dgrad[100][250])
BvalZ=antihelm(1.13, [0, 0, z])
gradZ=np.gradient(BvalZ[0], z)
Bhelm=helm(1.75, [0, 0, z])
#print max(BvalZ)
#Cgrad=appgrad(0.5, d)
linB=linearField(1.13, [0, 0, z])
#for i in R:
 #   temp=antihelm(i, [0, 0, z])
  #  Bz[i]=temp[0]
   # GradZ[i]=np.gradient(Bz[i])
    #DGradZ[i]=np.gradient(GradZ[i])
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('zero')
#ax.spines['right'].set_color('none')
#ax.spines['top'].set_color('none')
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')

# plot the function

#plt.plot(z, gradZ, 'r')
#plt.plot(z, Bhelm[0], 'b')

#plt.plot(z, np.gradient(np.gradient(BvalZ[0])), 'r')
#plt.plot(R, DGradZ, 'r')

#plt.plot(z,np.gradient(BvalZ[0]), 'r', z, np.gradient(BvalZ[1]), 'b')
Z=np.linspace(-0.04,0.04,50000)
Bfin=np.zeros(50000)
len=0.0
for i in range(0,10):
    d=0.002*i+0.020
    for j in range(0,10):
        rad=0.040-0.002*j
        a=rad/d
        z2=Z/d
        len=len+2*np.pi*rad
        Bfin=Bfin+adjustB2(a, z2, d)
print(len)
z = np.linspace(-1,1,1000)
Bnew=adjustB2(0.71, z,  0.045)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#for i in range(0,1500):
#    a0=i/500+0.25
#    Bplot=adjustB(a0, z, 0.05)
#    plt.plot(z, np.gradient(Bplot,z), 'r')
# plot the function
plt.plot(Z, Bfin, 'r')
#plt.plot(z, Bnew, 'r')
#plt.plot(Z, np.gradient(Bfin, Z), 'r')
print (Bfin[5001])
"""
# show the plot
plt.show()