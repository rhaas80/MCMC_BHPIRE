import corner
import numpy as np
import matplotlib.pyplot as plt      

chains=np.genfromtxt('chains.dat')

fig1 = plt.clf()
fig1 = corner.corner(chains, labels=[r"F$_1$",r"$\sigma_1$",r"x$_2$",r"y$_2$",r"F$_2$",r"$\sigma_2$"],
                     show_titles=True,title_fmt='.3f', bins=25, color='blue')

fig1.savefig('cornerplot.pdf')
