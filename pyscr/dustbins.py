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
ymax=5.0

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

subplot(3,4,1)
pln1=np.loadtxt('planetesimal.dat',skiprows=1)
pln=np.log10(pln1.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, pln, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title(' Planetesimal' ,fontsize=10)

sigd1=np.loadtxt('sigmadust.dat',skiprows=1)

subplot(3,4,2)
sigd=sigd1[:,0]+sigd1[:,1]+sigd1[:,2]+sigd1[:,3]+sigd1[:,4]+ \
   sigd1[:,5]+sigd1[:,6]+sigd1[:,7]+sigd1[:,8]+sigd1[:,9]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
#plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title(' Dust Total ',fontsize=10)

subplot(3,4,3)
sigd=sigd1[:,0]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('50 A Dust ',fontsize=10)

subplot(3,4,4)
sigd=sigd1[:,1]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('0.025um Dust ',fontsize=10)
plt.colorbar()

subplot(3,4,5)
sigd=sigd1[:,2]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('0.125um Dust ',fontsize=10)


subplot(3,4,6)
sigd=sigd1[:,3]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('0.6um Dust ',fontsize=10)

subplot(3,4,7)
sigd=sigd1[:,4]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('3.1um Dust ',fontsize=10)

subplot(3,4,8)
sigd=sigd1[:,5]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('15.8um Dust ',fontsize=10)
plt.colorbar()

subplot(3,4,9)
sigd=sigd1[:,6]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('79.4um Dust ',fontsize=10)

subplot(3,4,10)
sigd=sigd1[:,7]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('397um Dust ',fontsize=10)

subplot(3,4,11)
sigd=sigd1[:,8]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('2mm Dust ',fontsize=10)

subplot(3,4,12)
sigd=sigd1[:,9]
sigd=np.log10(sigd.reshape(nt,nr)+tiny)
levels = arange(-4,1,0.1)
plt.contourf(r, t, sigd, levels,cmap=plt.cm.gnuplot,extend='both')
plt.xlabel('log R ')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.ylabel('Time')
plt.title('1cm Dust ',fontsize=10)
plt.colorbar()


plt.savefig("dust.png")
#plt.show()



