import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from pylab import *
from matplotlib import animation


w1=0
w2=3
v1=0
v2=1

def grid(x, y, z, resX=100, resY=100):
#        xi=linspace(min(x),max(x),resX)
        xi=linspace(w1,w2,resX)
        yi=linspace(v1,v2,resY)
#        yi=linspace(min(y),max(y),resY)
        ZF=griddata(x,y,z,xi,yi,masked=True)
        X, Y = meshgrid(xi, yi)
        return X, Y, ZF

fig=plt.figure()
ax = fig.add_subplot(111,autoscale_on=False,xlim=(w1,w2),ylim=(v1,v2))

#time_text=ax.text(1,2.5,'',transform=ax.transAxes)
time_text=ax.text(1,2.5,'')

def init():
        time_text.set_text('')
        plt.xlabel('Radius')
        plt.ylabel('Height')
        plt.title(' Density (log$_{10}$ ; cm$^{-3}$)' )
        return[]

# Animation
nFrames=2
def animate(i):
# Disk parameters
        dfile="diskstruct_"+str(i*10+1)+".dat"
        disk = np.loadtxt(dfile,skiprows=2)
        r=disk[:,0]
        z=disk[:,1]/r
#        r = np.log10(disk[:,0]/1.5e13)
#        z = np.log10(disk[:,1]/1.5e13)
        r=np.log10(r/1.5e13)
        den = disk[:,2]/2.3e-24
        den = np.log10(den)
# Figure stuff
        levels = arange(3.0,13.0,1.0)
        X, Y, ZF = grid(r,z,den)
        cont=plt.contourf(X, Y, ZF, levels, cmap=plt.cm.hot,extend='both')
        time_text.set_text("t="+str(i+1))
        if i<1:
           plt.colorbar()
        return cont, time_text


anim=animation.FuncAnimation(fig,animate,interval=1,frames=nFrames,repeat=False,
         init_func=init)
#anim.save('den.mp4')
plt.show()


