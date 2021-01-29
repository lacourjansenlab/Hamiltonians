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

def visHam(H,N):
    # Clear Diagonal
    for i in range(N):
        H[i,i]=0.0
    plt.imshow(H)
    plt.axis('equal')
    plt.xlabel('Site')
    plt.ylabel('Site')
    plt.show()
    print(np.max(H))
    print(np.min(H))
