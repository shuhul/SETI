import numpy as np
import matplotlib.pyplot as plt
from pylab import *

f=np.loadtxt('junk.dat')
plt.plot(np.log10(f[:,0]),np.log10(f[:,1]))
plt.xlim(-3,4)
plt.ylim(-5,1)
plt.savefig('junk.png')

