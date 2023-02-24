import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time

from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
from pylab import *

###################################
#
# CONFIGURE
#
###################################
plt.rcParams['font.size'] = 8

###################################
#
# FUNCTIONS
#
###################################
def show_image(img,heading,vmin=2,vmax=3):
    cmap = plt.cm.gnuplot
    extent = [r[0],r[-1],t[0],t[-1]]

    plt.imshow(img,cmap=cmap,extent=extent,aspect=0.18,vmin=vmin,vmax=vmax,origin="lower")
    plt.xlim(0,3.0)
    plt.ylim(0,9.5)
    plt.ylabel('Time')
    plt.title(heading,fontsize=10)
    plt.colorbar()

###################################
#
# MAIN PROGRAM
#
###################################
tiny = 1.0e-30

rgrid = pd.read_csv('rgrid.inp').values[:,0]
r = np.log10(rgrid/1.5e13)
nr = r.size

tgrid = pd.read_csv('time.dat', names=["t"]).values[:,0]
t = tgrid/3.15e13
nt = t.size

plt.subplots_adjust(hspace=0.4)
plt.suptitle('Evolution')

subplot(3,2,1)

sigma1 = pd.read_csv('sigma.dat').values[:,0]
sigma = np.log10(sigma1.reshape(nt,nr)+tiny)
show_image(sigma,'Gas Surface Density',vmin=-2,vmax=3)

subplot(3,2,2)
g2d1 = pd.read_csv('gastodust.dat').values[:,0]
g2d = g2d1.reshape(nt,nr)
show_image(g2d,'Gas/Dust Ratio',vmin=0,vmax=100)

subplot(3,2,3)
pln1 = pd.read_csv('planetesimal.dat').values[:,0]
pln = np.log10(pln1.reshape(nt,nr)+tiny)
show_image(pln,'Planetesimal Density',vmin=-4,vmax=1)

##############
#
# MASS BINS.    sum(axis=0) == Sum rows.
#               sum(axis=1) == Sum columns.
#
subplot(3,2,4)

sigd1 = pd.read_csv('sigmadust.dat', skiprows=1, delim_whitespace=True, names=range(10))
sigd1 = sigd1.values
sigd = np.log10(sigd1.sum(axis=1).reshape(nt,nr)+tiny)

show_image(sigd,'Dust Surface Density',vmin=-4,vmax=1)

#
# Same `sigd1` as before.
#
subplot(3,2,5)
sigd = np.log10(sigd1[:,4].reshape(nt,nr)+tiny)

show_image(sigd,'3um Dust Surface Density',vmin=-4,vmax=1)
plt.xlabel('log R ')

#
# Same `sigd1` as before.
#
subplot(3,2,6)
sigd = np.log10(sigd1[:,8].reshape(nt,nr)+tiny)

show_image(sigd,'2mm Dust Surface Density',vmin=-4,vmax=1)
plt.xlabel('log R ')

plt.savefig("sigma-dc.png")
