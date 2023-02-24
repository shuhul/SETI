import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import sys 

timeframes=int(sys.argv[1])
deltat=int(sys.argv[2])

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

# Read alpha
dr=np.loadtxt('alpha.dat',skiprows=1)
alp=dr.reshape(nt,nr)

# Read planetesimals
dr=np.loadtxt('planetesimal.dat',skiprows=1)
pln=dr.reshape(nt,nr)

# Read Ices
dr=np.loadtxt('ices.dat',skiprows=1)
ices=dr[:,0].reshape(nt,nr)
vapor=dr[:,1].reshape(nt,nr)

# Read gas to dust
dr=np.loadtxt('gastodust.dat',skiprows=1)
gtd=dr.reshape(nt,nr)

# Read Dust Sigma 
dr=np.loadtxt('sigmadust.dat',skiprows=1)
s1=dr[:,1]
s6=dr[:,6]
s9=dr[:,9]
stot=dr[:,0]+s1+dr[:,2]+dr[:,3]+dr[:,4]+dr[:,5]+dr[:,6]+dr[:,7]+dr[:,8]+s9
s1=s1.reshape(nt,nr)
s6=s6.reshape(nt,nr)
s9=s9.reshape(nt,nr)
stot=stot.reshape(nt,nr)

s1r=np.zeros(nr)
s6r=np.zeros(nr)
s9r=np.zeros(nr)
sr=np.zeros(nr)
rad=np.zeros(nr)
spln=np.zeros(nr)
sd=np.zeros(nr)
sgtd=np.zeros(nr)
sice=np.zeros(nr)
svap=np.zeros(nr)

fig=plt.figure()
ax=plt.axes(xlim=(-1.5,4.0),ylim=(-10,7))
line_1, =  ax.plot([],[],c='teal')
line_2, =  ax.plot([],[],c='lime')
line_3, =  ax.plot([],[],c='r')
line_4, = ax.plot([],[],lw=3,c='black')
line_5, = ax.plot([],[],lw=3,c='green')
line_6, =  ax.plot([],[],c='magenta')
line_7, =  ax.plot([],[],ls='--',c='blue')
line_8, =  ax.plot([],[],ls='--',c='aqua')
time_text=ax.text(2.2,6.5,'')
label1_text=ax.text(2.2,6.0,' sub-um dust ',color='teal')
label2_text=ax.text(2.2,5.5,' micron dust',color='lime')
label3_text=ax.text(2.2,5.0,' mm dust',color='r')
label4_text=ax.text(2.2,4.5,' Gas    ',color='black')
label5_text=ax.text(2.2,4.0,' Pln    ',color='green')
label6_text=ax.text(2.2,3.5,' All dust',color='magenta')
label7_text=ax.text(2.2,3.0,' Water ice',color='blue')
label8_text=ax.text(2.2,2.5,' Water vapor',color='aqua')

def init():
   line_1.set_data([],[])
   line_2.set_data([],[])
   line_3.set_data([],[])
   line_4.set_data([],[])
   line_5.set_data([],[])
   line_6.set_data([],[])
   line_7.set_data([],[])
   line_8.set_data([],[])
   time_text.set_text('')
   return line_1,line_2,line_3,line_4,line_5,line_6,line_7,line_8,time_text

# 250(1), 409(10), 567(100), 643(300)
def stime(i):
    for j in range(0,1100):
         s1r[j]=s1[i,j]+tiny
         s6r[j]=s6[i,j]+tiny
         s9r[j]=s9[i,j]+tiny
         sr[j]=sall[i,j]+tiny
         spln[j]=pln[i,j]+tiny
         sd[j]=stot[i,j]+tiny
         rad[j]=r[j]+tiny
         sgtd[j]=gtd[i,j]+tiny
         sice[j]=ices[i,j]+tiny
         svap[j]=vapor[i,j]+tiny
         line_1.set_data(rad,np.log10(s1r))  
         line_2.set_data(rad,np.log10(s6r))  
         line_3.set_data(rad,np.log10(s9r))
         line_4.set_data(rad,np.log10(sr))
         line_5.set_data(rad,np.log10(spln))
         line_6.set_data(rad,np.log10(sd))
         line_7.set_data(rad,np.log10(sice))
         line_8.set_data(rad,np.log10(svap))
    time_text.set_text('t=%.2f Myr'% t[i])
    return line_1,line_2,line_3,line_4,line_5,line_6,line_7,line_8,time_text

anim=animation.FuncAnimation(fig,stime,init_func=init,blit=False,frames=timeframes,interval=deltat)

plt.show()



