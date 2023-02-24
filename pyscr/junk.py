import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

tiny=1.0e-30

#  Read R grid
rgrid=np.loadtxt('rgrid.inp',skiprows=1)
r=np.log10(rgrid/1.5e13)
nr=r.size

# Read times
tgrid=np.loadtxt('time.dat',skiprows=1)
t=tgrid/3.15e13
nt=t.size

# Read Gas Sigma
dr=np.loadtxt('sigma.dat',skiprows=1)
sall=dr.reshape(nt,nr)

# Read Gastodust
dr=np.loadtxt('gastodust.dat',skiprows=1)
g2d=dr.reshape(nt,nr)

# Read Gas Temperature
dr=np.loadtxt('temperature.dat',skiprows=1)
cs2=1.38e-16*dr/3.8e-24
cs2=cs2.reshape(nt,nr)

# Set dust sizes [in microns]
ad=[0.005,0.025,0.12,0.63,3.158,15.832,79.37,397.9,1994.7,10000.0]

# Read Dust Sigma and compute Stokes
#dr=np.loadtxt('sigmadust.dat',skiprows=1)
s1=np.zeros(nr)
s6=np.zeros(nr)
s9=np.zeros(nr)
stdr=np.zeros(nr)
stfrag=np.zeros(nr)
z=np.zeros(nr)
zs=np.zeros(nr)

# Make figure
fig=plt.figure()
ax=plt.axes(xlim=(-3.5,1.5),ylim=(-2,0))
scat, =ax.scatter([],[],c='r')
#line_1, =  ax.plot([],[],lw=2,c='cyan')
#line_2, =  ax.plot([],[],lw=2,c='black')
#line_3, =  ax.plot([],[],lw=2,c='r')
line_4, = ax.plot([],[],lw=3,c='blue')
line_5, = ax.plot([],[],lw=2,c='g')
line_6, = ax.plot([],[],lw=2,ls='--',c='g')
time_text=ax.text(0.8,-0.1,'')
label1_text=ax.text(0.5,-1.7,' sub-um dust ',color='cyan')
label2_text=ax.text(0.5,-1.8,' micron dust',color='black')
label3_text=ax.text(0.5,-1.9,' mm dust',color='r')

def init():
#   line_1.set_data([],[])
#   line_2.set_data([],[])
#   line_3.set_data([],[])
   xa=np.linspace(-3.5,1.5,50)
   za=-1.86+0.3*(xa+0.98)**2
   line_4.set_data(xa,za)
   time_text.set_text('')
   return line_4,time_text


def stime(i):
    for j in range(0,nr-1):
#         s1[j]=np.pi*ad[0]*1.0e-4*3.0/(2.0*p)
#         s6[j]=np.pi*ad[5]*1.0e-4*3.0/(2.0*p)
         z[j]=np.log10(1.0/(1.0+g2d[i,j]))
         p=sall[i,j]+0.0001
         s1[j]=np.pi*ad[0]*1.0e-4*3.0/(2.0*p)
         s6[j]=np.pi*ad[5]*1.0e-4*3.0/(2.0*p)
         s9[j]=np.pi*ad[8]*1.0e-4*3.0/(2.0*p)
         p=cs2[i,j]
#         stdr[j]=1.0e3*1.149e13/(np.sqrt(rgrid[j])*3.0*0.5*p)
         stdr[j]=0.1*1.0e3/300
         et_data(np.log10(s1),z)  
#         line_2.set_data(np.log10(s6),z)  
#         line_3.set_data(np.log10(s9),z)
         scat.set_array(np.log10(s1),z)
         line_5.set_data(np.log10(stdr),z)
         zs[j]=np.sqrt(1.0e-4/(1.0e-4+s9[j]))
         line_6.set_data(np.log10(zs),np.log10(s9))
    time_text.set_text('t=%.2f Myr'% t[i])
    return scat,line_5,line_6,time_text

anim=animation.FuncAnimation(fig,stime,init_func=init,blit=False,frames=400,interval=20)

plt.show()
