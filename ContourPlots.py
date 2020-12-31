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
    # B[0] and B[1] give Axial and Radial Magnetic Fields respectively in inverse units of
    #                       (mu I)/(2 pi d)

    x=r[0]
    y=r[1]
    z=r[2]
    rho=np.sqrt(x**2+y**2)    #radial distance in units of d
    B=[0, 0]
    k=np.sqrt((4*a*rho)/((a+rho)**2+(z-1)**2))
    K=sp.ellipk(k**2)
    E=sp.ellipe(k**2)
    B[0]=(1/np.sqrt((a+rho)**2+(z-1)**2))*(K+E*(a**2-rho**2-(z-1)**2)/((a-rho)**2+(z-1)**2))
    #if rho.all()==0:
        #B[1]=0
    #else:
    B[1]=((z-1)/rho)*(1/np.sqrt((a+rho)**2+(z-1)**2))*(-K+E*(a**2+rho**2+(z-1)**2)/((a-rho)**2+(z-1)**2))
    return B

def antihelm(a, r):
    "Anti-Helmholtz Coils"
    # The net magnetic field is difference of fields of two single coils
    # The single coils are kept at positions +d and -d from center
    # Equivalently, we can say that the axial coordinates of the single coil
    # are +z and -z in units of d, d>0 
    Bres=[0, 0]
    x=r[0]
    y=r[1]
    z=r[2]
    Bres1=singlecoil(a, [x, y, z])
    Bres2=singlecoil(a, [x, y, -z])
    Bres[0]=Bres1[0]-Bres2[0]
    Bres[1]=Bres1[1]+Bres2[1]
    return Bres

def adjustB_ax(a, x, z, d):
    "Ajusted Field"
    # For Anti-Helmholtz Coils, primary requirement is to return Gradient
    # Axial Magnetic Field variation along its Axis is found
    # The values of d are provided
    # The Magnetic Field is adjusted to return Gradient in the
    # units of Gauss/cm/A
    Badj=2*10**(-5)*antihelm(a, [x, 0, z])[0]/d
    return Badj

def adjustB_rad(a, x, z, d):
    "Ajusted Field"
    # For Anti-Helmholtz Coils, primary requirement is to return Gradient
    # Axial Magnetic Field variation along its Axis is found
    # The values of d are provided
    # The Magnetic Field is adjusted to return Gradient in the
    # units of Gauss/cm/A
    Badj2=2*10**(-5)*antihelm(a, [x, 0, z])[1]/d
    return Badj2

#Call the function

size=50
Z=np.linspace(-0.05,0.05,size) # Our Axis in units of meters
X=np.linspace(-0.05,0.05,size)
Bfin=np.zeros((size, size))
Bfin1=np.zeros((size, size))
Bfin2=np.zeros((size, size))
len=0.0  # To calculate the length of wire required in meters  

# Thickness of Copper wires is taken to be 2mm
# Each turn is taken as a separate Anti-Helmholtz Coil
# Summing over each turn gives us the total Magnetic Field Variation
# This method ensures that the calculation error is minimized
# Values of d and rad are taken for the middle of the copper wire of each turn
# i and j represent axial and radial turns respectively

for m in range(size):
    xval=X[m]
    #print(xval)
    for i in range(0,14):
        # Variation of distance of coil from the center
        d=0.025+0.002*i        
        for j in range(0,11):
            # Variation of radii of coils
            rad=0.025+0.002*j
            # Conversion into units of d
            a=rad/d
            z2=Z/d
            x2=xval/d
            len=len+2*np.pi*rad
            Bfin1[m]=Bfin1[m]+adjustB_ax(a, x2, z2, d)
            #print(Bfin1[m])
            Bfin2[m]=Bfin2[m]+adjustB_rad(a, x2, z2, d)
            Bfin[m]=np.sqrt(Bfin1[m]**2+Bfin2[m]**2)
print(len)  # Returns length of copper wire needed
# Plotting the Gradient
#print(Bfin)
#fig, ax = plt.subplots()
#CS = ax.contour(X, Z, Bfin)
#ax.clabel(CS, inline=1, fontsize=10)
#ax.set_title('Simplest default with labels')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#plt.contour([X, Z], Bfin)
#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('zero')
#ax.spines['right'].set_color('none')
#ax.spines['top'].set_color('none')
#ax.xaxis.set_ticks_position('bottom')
#ax.xaxis.set_label_position('bottom')
#ax.yaxis.set_ticks_position('left')
ax.set_xlabel('Axial Distance from center of system (meters)')
ax.set_ylabel('Radial Distance from center of system (meters)')
#ax.set_xlim(-0.0035, 0.0035)
#ax.set_ylim(7.0, 8.0)
#ax.set_ylim(0.0, 10.0)
plt.contourf(X, Z, Bfin, cmap='gray')
cbar = plt.colorbar()
#cbar.ax.set_yticklabels(['0','1','2','>3'])
cbar.set_label('Magnetic Field in Gauss/A')


# As a final step, we plot gradient as a function of 
# Axial Distance from center in meters
#plt.plot(Z, Bfin1[0], 'k')

# show the plot
plt.show()