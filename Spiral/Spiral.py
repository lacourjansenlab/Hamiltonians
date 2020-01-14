import numpy as np
import matplotlib.pyplot as plt

# Global parameters
d2eA=0.20819434 # Debye to eÅ
bohr3=0.529177210 # bohr in Å
Eh2icm=219500 # Hartree to cm-1
A=Eh2icm*bohr3*d2eA**2
print(A)

# Define Parameters
d=16.0 # Distance
r=3.0 # Radius
N=300 # Number of molecules
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

# Plot the structure
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
    J=np.inner(mu[n,:],mu[m,:])/d3+np.inner(mu[n,:],dx)*np.inner(dx,mu[m,:])/d5
    H[n,m]=J*A
    H[m,n]=J*A
#print(H)
E,c=np.linalg.eig(H)
print(E)

# Make spectrum
bins=1000
EminA=np.min(E)
EmaxA=np.max(E)
Emin=EminA-0.1*(EmaxA-EminA)
Emax=EmaxA+0.1*(EmaxA-EminA)
dE=(Emax-Emin)/bins
Ex=np.linspace(Emin,Emax,bins)
Ey=np.zeros(bins)
for n in range(N):
  bin=int(round((E[n]-Emin)/dE))
  Emu=np.zeros(3)
  Emu[0]=np.inner(c[:,n],mu[:,0])
  Emu[1]=np.inner(c[:,n],mu[:,1])
  Emu[2]=np.inner(c[:,n],mu[:,2])
  Ey[bin]=Ey[bin]+np.linalg.norm(Emu)**2

plt.plot(Ex,Ey)
plt.show()

# Convolute spectrum
# First create normalized Gaussian centered in the middle
sigma=10
Cx=Ex-(Emax+Emin)/2 # Make new axis with value zero in the middle of array
Cy=np.exp(-Cx**2/2/sigma**2)/np.sqrt(2*np.pi*sigma**2)
plt.plot(Cx,Cy)
plt.show()

# Do the actual convolusion
Ny=np.convolve(Ey,Cy,mode='same')
plt.plot(Ex,Ny)
plt.show()

# Plot everything in one final plot
plt.plot(Ex,Ey/np.max(Ey))
plt.plot(Ex,Cy/np.max(Cy))
plt.plot(Ex,Ny/np.max(Ny))
plt.show()
