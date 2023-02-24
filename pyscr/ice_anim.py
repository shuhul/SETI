import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import click
import util


default_timeframes = int(util.load('time.info', 0))

timeframes = click.prompt("Timeframes", type=int, default=default_timeframes)
deltat = click.prompt("Deltat", type=int, default=1)

tiny = 1.0e-30

#  Read R grid
rgrid = util.load('rgrid.inp', 1)
r = np.log10(rgrid/1.5e13)
# r=(rgrid/1.5e13)
nr = r.size

# Read times
tgrid = util.load('time.dat', 1)
t = tgrid/3.15e13
nt = t.size

# Read Gas Sigma
dr = util.load('sigma.dat', 2)
sall = dr.reshape(nt, nr)


# Read planetesimals
dr = util.load('planetesimal.dat', 2)
pln = dr.reshape(nt, nr)

# Read Ices
dr = util.load('ices.dat', 2)
# dr = dr *1e3
wice = dr[:, 0].reshape(nt, nr)
wvapor = dr[:, 1].reshape(nt, nr)
coice = dr[:, 2].reshape(nt, nr)
covapor = dr[:, 3].reshape(nt, nr)

# Read gas to dust
dr = util.load('gastodust.dat', 2)
gtd = dr.reshape(nt, nr)

# Read Dust Sigma
dr = util.load('sigmadust.dat', 2)
s1 = dr[:, 1]
s6 = dr[:, 6]
s9 = dr[:, 9]
stot = dr[:, 0]+s1+dr[:, 2]+dr[:, 3]+dr[:, 4] + \
    dr[:, 5]+dr[:, 6]+dr[:, 7]+dr[:, 8]+s9
s1 = s1.reshape(nt, nr)
s6 = s6.reshape(nt, nr)
s9 = s9.reshape(nt, nr)
stot = stot.reshape(nt, nr)

rad = np.zeros(nr)
spln = np.zeros(nr)
sd = np.zeros(nr)
sr = np.zeros(nr)
sgtd = np.zeros(nr)
swice = np.zeros(nr)
swvap = np.zeros(nr)
scoice = np.zeros(nr)
scovap = np.zeros(nr)
smmdust = np.zeros(nr)

fig = plt.figure()
ax = plt.axes(xlim=(-0.5, 2.5), ylim=(-8, 5))
line_1, =  ax.plot([], [], lw=3, c='k')
line_2, =  ax.plot([], [], lw=3, c='g')
line_3, =  ax.plot([], [], lw=3, c='magenta')
line_4, = ax.plot([], [], c='blue')
line_5, = ax.plot([], [], c='aqua', ls=':', lw=2)
line_6, =  ax.plot([], [], c='red')
line_7, =  ax.plot([], [], ls=':', c='magenta', lw=2)
# line_8, =  ax.plot([],[],c='red')
time_text = ax.text(2.2, 6.5, '')
label1_text = ax.text(2.2, 4.5, ' Gas    ', color='black')
label2_text = ax.text(2.2, 4.0, ' Pln    ', color='green')
label3_text = ax.text(2.2, 3.5, ' Dust', color='magenta')
label4_text = ax.text(2.2, 3.0, ' Water ice', color='blue')
label5_text = ax.text(2.2, 2.5, ' Water vapor', color='aqua')
label6_text = ax.text(2.2, 2.0, ' CO ice', color='red')
label7_text = ax.text(2.2, 1.5, ' CO gas', color='magenta')
# label8_text=ax.text(2.2,1.0,' mmdust',color='red')


def init():
    line_1.set_data([], [])
    line_2.set_data([], [])
    line_3.set_data([], [])
    line_4.set_data([], [])
    line_5.set_data([], [])
    line_6.set_data([], [])
    line_7.set_data([], [])
#   line_8.set_data([],[])
    time_text.set_text('')
    return line_1, line_2, line_3, line_4, line_5, line_6, line_7, time_text

# 250(1), 409(10), 567(100), 643(300)


def stime(i):
    for j in range(0, 1100):
        sr[j] = sall[i, j]+tiny
        spln[j] = pln[i, j]+tiny
        sd[j] = stot[i, j]+tiny
        rad[j] = r[j]+tiny
        scoice[j] = coice[i, j]+tiny
        scovap[j] = covapor[i, j]+tiny
        swice[j] = wice[i, j]+tiny
        swvap[j] = wvapor[i, j]+tiny
        smmdust[j] = s9[i, j] + tiny
        # scoice[j]=(coice[i,j]+tiny)/(coice[i,j]+wice[i,j]+tiny)
        # scovap[j]=(covapor[i,j]+tiny)/(covapor[i,j]+wvapor[i,j]+tiny)
        line_1.set_data(rad, util.log10(sr))
        line_2.set_data(rad, util.log10(spln))
        line_3.set_data(rad, util.log10(sd))
        line_4.set_data(rad, util.log10(swice))
        line_5.set_data(rad, util.log10(swvap))
        line_6.set_data(rad, util.log10(scoice))
        line_7.set_data(rad, util.log10(scovap))
#         line_8.set_data(rad,np.log10(smmdust))
    time_text.set_text('t=%.2f Myr' % t[i])
    return line_1, line_2, line_3, line_4, line_5, line_6, line_7, time_text


anim = animation.FuncAnimation(
    fig, stime, init_func=init, blit=False, frames=timeframes, interval=deltat)

plt.show()
