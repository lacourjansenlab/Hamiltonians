import numpy as np
import matplotlib.pyplot as plt

### Define Plotting Defaults
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

# This subroutine calculates stick spectra and spectra convolouted with
# a Gaussian. Provide the Hamiltonian, H, the transition-dipoles, mu,
# the size of the Hamiltonian, N, and the standard deviation for the
# Gaussian convolution, sigma. Plots are generated. This subroutine is
# intended as a way to test the Hamiltonians generated in this package.
def CD(H,mu,r,N,sigma):
  # Diagonalize Hamiltonian
  E,c=np.linalg.eigh(H)

  # Make spectrum
  bins=1000
  EminA=np.min(E)
  EmaxA=np.max(E)
  # Ensure that the resolution is bigger than the witdh of each peak
  if (EmaxA-EminA)/sigma>1000:
      bins=int(np.floor((EmaxA-EminA)/sigma))
  # Add 10% of full bandwidth and three times the convolution width on each side
  Emin=EminA-0.1*(EmaxA-EminA)-3*sigma
  Emax=EmaxA+0.1*(EmaxA-EminA)+3*sigma
  dE=(Emax-Emin)/bins
  Ex=np.linspace(Emin,Emax,bins)
  Ey=np.zeros(bins)
  for n in range(N):
    bin=int(round((E[n]-Emin)/dE))
    Emu=np.zeros(3)
    R=0.0
    for k in range(N):
        for m in range(N):
            R=R+c[k,n]*c[m,n]*mu[k,0]*mu[m,1]*(r[k,2]-r[m,2])
            R=R+c[k,n]*c[m,n]*mu[k,1]*mu[m,2]*(r[k,0]-r[m,0])
            R=R+c[k,n]*c[m,n]*mu[k,2]*mu[m,0]*(r[k,1]-r[m,1])
            R=R-c[k,n]*c[m,n]*mu[k,0]*mu[m,2]*(r[k,1]-r[m,1])
            R=R-c[k,n]*c[m,n]*mu[k,2]*mu[m,1]*(r[k,0]-r[m,0])
            R=R-c[k,n]*c[m,n]*mu[k,1]*mu[m,0]*(r[k,2]-r[m,2])
    Ey[bin]=Ey[bin]+R

  plt.plot(Ex,Ey)
  plt.xlabel('Wavenumbers [cm$^{-1}$]',fontsize=16)
  plt.ylabel('CD [arb.u.]', fontsize=16)
  plt.show()

  # Convolute spectrum
  # First create normalized Gaussian centered in the middle
  Cx=Ex-(Emax+Emin)/2 # Make new axis with value zero in the middle of array
  Cy=np.exp(-Cx**2/2/sigma**2)/np.sqrt(2*np.pi*sigma**2)
  plt.plot(Cx,Cy)
  plt.xlabel('Wavenumbers [cm$^{-1}$]',fontsize=16)
  plt.ylabel('CD [arb.u.]', fontsize=16)
  plt.show()

  # Do the actual convolusion
  Ny=np.convolve(Ey,Cy,mode='same')
  plt.plot(Ex,Ny)
  plt.show()

  # Plot everything in one final plot
  plt.plot(Ex,Ey/np.max(Ey))
  plt.plot(Ex,Cy/np.max(Cy))
  plt.plot(Ex,Ny/np.max(Ny))
  plt.xlabel('Wavenumbers [cm$^{-1}$]',fontsize=16)
  plt.ylabel('CD [arb.u.]', fontsize=16)
  plt.show()
