import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from pylab import *

plt.rcParams['font.size']=10

tiny=1.0e-30
xmin=0
xmax=3.0
ymin=0.0
ymax=5.3

def grid(x, y, z, resX=100, resY=100):
#        xi=linspace(min(x),max(x),resX)
#        yi=linspace(min(y),max(y),resY)
        xi=linspace(xmin,xmax,resX)
        yi=linspace(ymin,ymax,resY)
        Z=griddata(x,y,z,xi,yi)
        X, Y = meshgrid(xi, yi)
        return X, Y, Z

rgrid=np.loadtxt('rgrid.inp',skiprows=1)
r=np.log10(rgrid/1.5e13)
nr=r.size

tgrid=np.loadtxt('time.dat',skiprows=1)
t=tgrid/3.15e13
nt=t.size

g2smpy=np.log10(3.15e7/1.99e33)

plt.subplots_adjust(hspace=0.4)
subplot(2,2,1)
acc1=np.loadtxt('mdot.dat',skiprows=1)
acc=np.log10(acc1.reshape(nt,nr)+tiny)+g2smpy
levels = arange(-15,-5,1)
plt.contourf(r, t, acc, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel('Log Radius (AU) ')
plt.ylabel('Time in Myr')
plt.title(' Log Mdotacc ',fontsize=16)
plt.colorbar()

subplot(2,2,2)
sdot=np.loadtxt('sigmadotevap.dat',skiprows=1)
mpe=sdot.reshape(nt,nr)
levels = arange(0,100,5)
plt.contourf(r, t, g2d, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel('Log Radius (AU) ')
plt.ylabel('Time in Myr')
plt.title(' Gas/Dust Ratio ',fontsize=16)
plt.colorbar()

plt.show()

