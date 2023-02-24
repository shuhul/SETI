import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from pylab import *
from matplotlib import animation

def grid(x, y, z, resX=100, resY=100):
#        xi=linspace(min(x),max(x),resX)
        xi=linspace(0.0,3.0,resX)
        yi=linspace(min(y),max(y),resY)
        Z=griddata(x,y,z,xi,yi)
        X, Y = meshgrid(xi, yi)
        return X, Y, Z

# Making the figure
plt.subplots_adjust(hspace=0.4)
# subplot(numrows,numcols,fignum)
subplot(2,2,1)

# Animation
def animate(i):
# Disk parameters
        dfile='diskstruct_3.dat'
        disk = np.loadtxt(dfile,skiprows=2)

        r = disk[:,0]
        z = disk[:,1]/r
        r=np.log10(r/1.5e13)
        den = disk[:,2]/2.3e-24
        den = np.log10(den)
        tem = np.log10(disk[:,3])
        levels = arange(0.0,4.25,0.25)
        X, Y, Z = grid(r,z,tem)
        contd=plt.contourf(X, Y, Z, levels, cmap=plt.cm.hot)
        plt.xlabel('Radius')
        plt.ylabel('Height')
        plt.title(' Temperature (log$_{10}$ ; K) ')
        plt.colorbar()

subplot(2,2,2)
levels = arange(3.0,13.0,1.0)
X, Y, Z = grid(r,z,den)
plt.contourf(X, Y, Z, levels, cmap=plt.cm.bone,extend='both')
plt.xlabel('Radius')
plt.ylabel('Height')
plt.title(' Density (log$_{10}$ ; cm$^{-3}$) ')
plt.colorbar()


        

# Show and save
#plt.savefig("disk.png")
plt.show()



