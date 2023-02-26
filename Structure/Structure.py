# This python script require vpython to be installed (see vpython.org)
import numpy as np
import matplotlib.pyplot as plt
from vpython import *

# This subroutine visualizes a 3D structure using vpython
# The vpython library must be available

# Visualize the bare structure with white arrows for the transition dipoles
def visual(x,mu,N,scale):
    for i in range(N):
      arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(scale*mu[i,0],scale*mu[i,1],scale*mu[i,2]),round=True)

# Visualize a selected exciton state with white arrows showing the phase and magnitude of the wave function. The molecules are shown as small spheres.
def visual_exciton(x,mu,c,index,N,scale):
    for i in range(N):
        # Show molecule
        sphere(pos=vector(x[i,0],x[i,1],x[i,2]),radius=0.1)
        # Show associated transition dipole
        arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(c[i,index]*scale*mu[i,0],c[i,index]*scale*mu[i,1],c[i,index]*scale*mu[i,2]),round=True)

# Visualize a selected exciton state with white arrows showing the phase and magnitude of the wave function. The molecules are shown as small spheres.
def visual_exciton_color(x,mu,c,index,N,scale):
    cmax=np.amax(np.abs(c[:,index]))
    for i in range(N):
        # Show molecule
        sphere(pos=vector(x[i,0],x[i,1],x[i,2]),radius=0.1)
        # Show associated transition dipole
        green=np.abs(1-c[i,index]*c[i,index]/(cmax**2))
        red=green
        blue=green 
        if c[i,index]<0:
          red=1
        if c[i,index]>0:
          blue=1
        col=vector(red,green,blue)
        print(i,col)
        arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(scale*mu[i,0],scale*mu[i,1],scale*mu[i,2]),color=col,round=True)

