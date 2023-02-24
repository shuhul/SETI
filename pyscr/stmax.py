import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from pylab import *

plt.rcParams['font.size']=8

tiny=1.0e-30
xmin=0.0
xmax=3.0
ymin=0.0
ymax=9.5

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

plt.subplots_adjust(hspace=0.4)
plt.suptitle('values less than alpha will drift')

ufrag=1000.0
sigg=np.loadtxt('sigma.dat',skiprows=1)
temp=np.loadtxt('temperature.dat',skiprows=1)
cs=sqrt(1.38e-16*temp/3.9e-24)
st=(3.5*np.pi*0.5*1.0/sigg)
# above is the Stokes number for a 1cm size grain
# Stmax = ufrag**2/(alpha*cs**2)
# assume amax=1cm, and this gives Stmax < 1 (no drift) implies
#    alpha > ufrag**2/cs**2 for no drift
alpha=(ufrag**2/cs**2)
alpha=np.log10(alpha.reshape(nt,nr)+tiny)
levels = arange(-6,-2,0.1)
plt.contourf(r, t, alpha, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title(' Crit alpha ',fontsize=10)
plt.colorbar()

plt.savefig("stmax.png")




