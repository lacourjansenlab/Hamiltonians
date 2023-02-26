import numpy as np
import matplotlib.pyplot as plt
# This python script require vpython to be installed (see vpython.org)
# vpython is used for the structural visualization
from vpython import *
import sys
sys.path.append("../Spectra")
sys.path.append("../Structure")
import Spectra
import CD
import Structure

# Global parameters
d2eA=0.20819434 # Debye to eÅ
bohr3=0.529177210 # bohr in Å
Eh2icm=219500 # Hartree to cm-1
A=Eh2icm*bohr3*d2eA**2
print(A)

# Define Parameters
r=20.0 # Radius of cylinder in Ångstrøm
N1=10 #Number of rings
N2=10 # Number of molecules in ring
N=N1*N2 # Total number of molecules
h=3 # Distance between rings in Ångstrøm
d2r=np.pi/180.0 # Degree to radians
alpha=4.0*d2r # Alpha angle for transition dipole
beta=55.0*d2r # Beta angle for transition dipole
delta=5*d2r # Turn for each ring
mum=5.5 # Dipole in Debye

# Create positions
x=np.zeros((N,3))
mu=np.zeros((N,3))
for n1 in range(N1):
  for n2 in range(N2):
    n=N2*n1+n2
    p=2*np.pi/N2*n2+n1*delta
    x[n,0]=r*np.cos(p)
    x[n,1]=r*np.sin(p)
    x[n,2]=h*n1
    mu[n,2]=mum*np.cos(beta)
    mu[n,0]=mum*np.sin(beta)*(-np.sin(p)*np.cos(alpha)+np.cos(p)*np.sin(alpha))
    mu[n,1]=mum*np.sin(beta)*(np.cos(p)*np.cos(alpha)+np.sin(p)*np.sin(alpha))
    lmu=np.linalg.norm(mu[n,:])
#  print(lmu)  

# Plot the structure in 2D
plt.plot(x[:,0],x[:,1])
for n in range(N):
  plt.arrow(x[n,0],x[n,1],mu[n,0],mu[n,1])
# Make x and y direction equivalent on screen
plt.axis('equal')
plt.show()

# Create Hamiltonian
H=np.zeros((N,N))
for n in range(N):
  for m in range(n+1,N):
    dx=x[n,:]-x[m,:]
    dd=np.linalg.norm(dx)
    d3=dd*dd*dd
    d5=d3*dd*dd
    J=np.inner(mu[n,:],mu[m,:])/d3-3*np.inner(mu[n,:],dx)*np.inner(dx,mu[m,:])/d5
    H[n,m]=J*A
    H[m,n]=J*A

# Plot structure
#Structure.visual(x,mu,N,1)
# Make spectrum
#Spectra.absorption(H,mu,N,10)
#CD.CD(H,mu,x,N,10)
# Visualize state
E,c=np.linalg.eigh(H)
Structure.visual_exciton_color(x,mu,c,1,N,1)

