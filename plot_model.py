import numpy as np                    # imports library for math
import matplotlib.pyplot as plt       # imports library for plots
import csv                            # import csv reader

fname='model.dat'

alldataRead = np.recfromcsv(fname, case_sensitive=True) # read all data
uCo=alldataRead['uCo']/1.e9
vCo=alldataRead['vCo']/1.e9
vis=np.absolute(alldataRead['VisAmp'])
error=alldataRead['Sigma']
model=alldataRead['Model']

b0=np.sqrt(uCo*uCo+vCo*vCo)            # baseline length
plt.errorbar(b0,vis,yerr=error,fmt='none',capsize=2,label='data')
plt.plot(b0,model,'ro',markersize=3,label='model')

plt.xlabel('Baseline length')
plt.ylabel('Visibility Amplitude')
plt.legend()

plt.yscale('log')
plt.show()
