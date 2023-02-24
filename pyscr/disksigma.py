import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from pylab import *

plt.rcParams['font.size']=12

tiny=1.0e-30
xmin=0.0
xmax=2.8
ymin=0.0
ymax=3.0

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

subplot(2,2,1)
sigma1=np.loadtxt('sigma.dat',skiprows=1)
sigma=np.log10(sigma1.reshape(nt,nr)+tiny)
levels = arange(-2,3,0.1)
plt.contourf(r, t, sigma, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('Log $\Sigma_{GAS}(r)$ ',fontsize=12)
plt.colorbar()

subplot(2,2,4)
g2d1=np.loadtxt('gastodust.dat',skiprows=1)
#g2d=np.log10(g2d1.reshape(nt,nr)+tiny)
g2d=g2d1.reshape(nt,nr)
#levels = arange(-2,2.1,0.1)
levels = arange(0,100,5)
plt.contourf(r, t, g2d, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel('Log R ')
plt.title(' Gas/Dust Ratio ',fontsize=12)
plt.colorbar()

subplot(2,2,3)
sigd1=np.loadtxt('sigmadust.dat',skiprows=1)
sigd=sigd1[:,0]+sigd1[:,1]+sigd1[:,2]+sigd1[:,3]+sigd1[:,4]+ \
   sigd1[:,5]+sigd1[:,6]+sigd1[:,7]+sigd1[:,8]+sigd1[:,9]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.xlabel('Log R ')
plt.title(' Log $\Sigma_{DUST}(r)$' ,fontsize=12)
plt.colorbar()


subplot(2,2,2)
#levels = arange(-2,2.1,0.1)
sdot=np.loadtxt('mdotpe.dat',skiprows=1)
sdot=np.log10(abs(sdot.reshape(nt,nr)+tiny))
levels = arange(-13,-6,0.1)
plt.contourf(r, t, sdot, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.title(' Mdot$_{PE}(r)$ ',fontsize=12)
plt.colorbar()

#subplot(3,2,6)
#sigd1=np.loadtxt('sigmadust.dat',skiprows=1)
#sigd=sigd1[:,8]
#sigd=np.log10(sigd.reshape(nt,nr)+tiny)
#levels = arange(-4,1,0.1)
#plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
#plt.xticks(fontsize=10)
#plt.yticks(fontsize=10)
#plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#plt.ylabel('Time')
#plt.title(' 2mm Dust Surface Density ',fontsize=10)
#plt.colorbar()

plt.savefig("sigmagasdust.png")
#plt.show()



