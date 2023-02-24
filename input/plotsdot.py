import numpy as np
import matplotlib.pyplot as plt


r=np.loadtxt('rgrid.inp',skiprows=1)
s=np.loadtxt('sigmadot.inp-pemodel',skiprows=1)
snew=np.loadtxt('sigmadot.inp',skiprows=1)


plt.loglog(r/1.5e13,s/1.0e-15,c='k')
plt.loglog(r/1.5e13,snew/1.0e-15,c='b')

plt.ylim(1.0e-3,1.0e3)


dr=np.diff(r)
delr = np.insert(dr,0,0)

mdot = np.sum(2.0*np.pi*delr*r*s)
mdotnew = np.sum(2.0*np.pi*delr*r*snew)
print(mdot,' orig ')
print(mdotnew,' new ')

plt.show()
