import numpy as np
from matplotlib import pyplot as plt
import sys


title=sys.argv[1]

fig,ax1=plt.subplots()

tiny=1.0e-30
d=np.loadtxt("diskmass.dat")
t=d[:,0]/1.0e6
mg=np.log10(d[:,1]+tiny)
md=np.log10(d[:,2]+tiny)
mp=np.log10(d[:,3]+tiny)
mwat=np.log10(d[:,4]+tiny)
mco=np.log10(d[:,5]+tiny)

ax1.plot(t,mg,'-',color='black')
ax1.plot(t,mg-2,'.',color='black')
ax1.plot(t,md,'-',color='r')
ax1.plot(t,mp,'-',color='g')
ax1.plot(t,mwat,'.',color='b')
ax1.plot(t,mco,'.',color='b')
ax1.set_xlim(0,9)
ax1.set_ylim(-10,0.1)
ax1.set_xlabel('Time in Myr')
ax1.set_ylabel(' Log Mass (M$_{\odot}$) ')
ax1.text(1,-1,' M$_{G}$(t) ',color='black')
ax1.text(1,-3,' M$_{D}$(t) ',color='red')
ax1.text(1,-4.5,' M$_{P}$(t) ',color='green')
ax1.text(4,-1,title,color='black',size=20)

#plt.savefig('mass.png')
#--------------------------------
#  Read R grid
rgrid=np.loadtxt('rgrid.inp',skiprows=1)
r=rgrid
nr=r.size
dr=np.zeros(nr)
dr[0]=2.0*np.pi*r[0]*r[0] # dr is r itself
for j in range(1,nr-1):
    dr[j]=np.pi*(r[j-1]+r[j])*(r[j]-r[j-1])
r=np.log10(r/1.5e13)   

# Read times
tgrid=np.loadtxt('time.dat',skiprows=1)
t=tgrid/3.15e13
nt=t.size

# Read planetesimals
pm=np.loadtxt('planetesimal.dat',skiprows=1)
pln=pm.reshape(nt,nr)
spln=np.zeros(nr)
for j in range(0,nr-1):
   dummy=pln[nt-1,j]*dr[j]/1.9e33 + tiny
   spln[j]=np.log10(dummy)
#--------------------------------

ax2=ax1.twiny()
ax2.plot(r,spln,'g--')
ax2.set_xlim(-1,3.3)
ax2.set_ylim(-10,0.1)
ax2.set_xlabel('Log Radius (AU)')
ax2.text(0.5,-9,' M$_{P}$(R) ',color='green')


plt.savefig('mass.png')





