import numpy as np
import matplotlib.pyplot as plt

# This subroutine calculates stick spectra and spectra convolouted with
# a Gaussian. Provide the Hamiltonian, H, the transition-dipoles, mu,
# the size of the Hamiltonian, N, and the standard deviation for the
# Gaussian convolution, sigma. Plots are generated. This subroutine is
# intended as a way to test the Hamiltonians generated in this package.
def absorption(H,mu,N,sigma):
  # Diagonalize Hamiltonian
  E,c=np.linalg.eig(H)

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
    Emu[0]=np.inner(c[:,n],mu[:,0])
    Emu[1]=np.inner(c[:,n],mu[:,1])
    Emu[2]=np.inner(c[:,n],mu[:,2])
    Ey[bin]=Ey[bin]+np.linalg.norm(Emu)**2

  plt.plot(Ex,Ey)
  plt.show()

  # Convolute spectrum
  # First create normalized Gaussian centered in the middle
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
