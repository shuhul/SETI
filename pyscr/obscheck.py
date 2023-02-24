import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import sys


fname=sys.argv[1]

d=np.loadtxt(fname)
time=d[:,0]
mdot=d[:,1]
mdust=d[:,2]
mgas=d[:,3]

m=np.loadtxt('manara16.txt')
mx=m[:,0]
my=m[:,1]

h=np.loadtxt('Hartmann98.txt')
hx=h[:,0]
hy=h[:,1]

plt.subplots_adjust(hspace=0.4)
plt.suptitle(fname)

subplot(2,2,1)
plt.scatter(hx,hy)
z=np.polyfit(hx,hy,1)
print z[1],z[0]
zy=z[0]*hx+z[1]
plt.scatter(hx,zy,c='r')
plt.plot(time,mdot)
plt.xlim(4.9,7.1)
plt.ylim(-12,-6)
plt.xlabel('Log Time (yr) ')
plt.ylabel('Log dM$_*$/dt (M$_{\odot}$/yr) ')


subplot(2,2,2)
plt.scatter(my,mx)
z=np.polyfit(my,mx,1)
print z[1],z[0]
zy=z[0]*my + z[1]
plt.scatter(my,zy,c='r')
x=mdot[(time > 6) & ( time < 6.5)]
y=mdust[(time > 6) & ( time < 6.5)]
z100=mgas[(time > 6) & ( time < 6.5)]
z50=z100+2.0
#plt.plot(mdot,mdust)
#plt.plot(mdot,mgas)
plt.plot(x,y,lw=3) 
plt.plot(x,z100,lw=3,c='g')
plt.plot(x,z50,lw=3,c='g')
plt.xlim(-12,-6)
plt.ylim(-7,-3)
plt.xlabel('Log dM$_*$/dt (M$_{\odot}$/yr) ')
plt.ylabel(' Log M$_{dust}$ (M$_{\odot}$) ')

plt.show()


