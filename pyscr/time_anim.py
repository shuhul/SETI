import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import sys 

timeframes=int(sys.argv[1])
fname=sys.argv[2]
valmin=float(sys.argv[3])
valmax=float(sys.argv[4])


tiny=1.0e-30

#  Read R grid
rgrid=np.loadtxt('rgrid.inp',skiprows=1)
r=np.log10(rgrid/1.5e13)
nr=r.size

# Read times
tgrid=np.loadtxt('time.dat',skiprows=1)
t=tgrid/3.15e13
nt=t.size

# Read file 
dr=np.loadtxt(fname,skiprows=1)
dr=dr[:,0]
sall=dr.reshape(nt,nr)

lpos=0.9*valmax
sr=np.zeros(nr)
rad=np.zeros(nr)

fig=plt.figure()
ax=plt.axes(xlim=(-0.5,4.0),ylim=(valmin,valmax))
line_4, = ax.plot([],[],lw=3,c='blue')
time_text=ax.text(2.2,lpos,'')
plt.title(fname,color='blue')

def init():
   line_4.set_data([],[])
   time_text.set_text('')
   return line_4,time_text

def stime(i):
    for j in range(0,1100):
         sr[j]=sall[i,j]+tiny
         rad[j]=r[j]+tiny
         line_4.set_data(rad,np.log10(abs(sr)))

    time_text.set_text('t=%.2f Myr'% t[i])
    return line_4,time_text

anim=animation.FuncAnimation(fig,stime,init_func=init,blit=False,frames=timeframes,interval=1)

plt.show()
