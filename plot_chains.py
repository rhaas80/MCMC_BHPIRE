import numpy as np                    # imports library for math
import matplotlib.pyplot as plt       # imports library for plots
import csv                            # import csv reader

samples=np.genfromtxt('chains.dat')

# for all parameters
for iparam in range(0,samples.shape[1]):
    yarray=samples[:,iparam]          # pick chain for one parameter
    yarray/=yarray[yarray.size-1]     # normalize to its last value
    xarray=np.arange(0,yarray.size)   # make a counting array
    plt.plot(xarray,yarray,label=str(iparam))  # plot them
    

plt.xlabel('Chain Number')
plt.ylabel('Normalized Parameter Value')
plt.legend()

#plt.yscale('log')
plt.show()
