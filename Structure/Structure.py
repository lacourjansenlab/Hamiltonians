# This python script require vpython to be installed (see vpython.org)
import numpy as np
import matplotlib.pyplot as plt
from vpython import *

# This subroutine visualizes a 3D structure using vpython
# The vpython library must be available

def visual(x,mu,N,scale):
    for i in range(N):
      arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(scale*mu[i,0],scale*mu[i,1],scale*mu[i,2]))

