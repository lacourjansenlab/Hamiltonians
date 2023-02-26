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
d=16.0 # Distance between molecules along spiral in Ångstrøm
r=3.0 # Distance between layers in Ångstrøm
N=20 # Number of molecules
n0=2 # Number of first point
d2r=np.pi/180.0 # Degree to radians
alpha=4.0*d2r # Alpha angle for transition dipole
beta=55.0*d2r # Beta angle for transition dipole
mum=5.5 # Dipole in Debye

# Create positions
x=np.zeros((N,3))
mu=np.zeros((N,3))
for n in range(N):
  p=np.sqrt(2*d*(n+n0)/r) # Leave our first n0 points
  x[n,0]=r*p*np.cos(p)
  x[n,1]=r*p*np.sin(p)
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
Structure.visual(x,mu,N,1)
# Make spectrum
Spectra.absorption(H,mu,N,10)
CD.CD(H,mu,x,N,10)
# Visualize state
#visual_exciton(x,mu,c,index,N,scale)

