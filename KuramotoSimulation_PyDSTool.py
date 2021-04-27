from PyDSTool import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import numpy.random as random

kuramotoDS=args()
##################################################################
###### Specify number of oscillators and strength of coupling here
##################################################################

num_oscillators=15
coupling_strength=0.009

kuramotoDS.pars={'omega':0.1,'K':coupling_strength}

kuramotoDS.varspecs={'x[i]':'for(i,0,'+str(num_oscillators-2)+',omega+K*sum(j,0,'+str(num_oscillators-2)+',sin(x[j]-x[i])))', 'x'+str(num_oscillators-1):'omega+K*sum(i,0,'+str(num_oscillators-1)+',sin(x[i]-x'+str(num_oscillators-1)+'))'}

kuramotoDS.tdata=[0,50*np.pi]
kuramotoDS.ics={}
varnamedict=[]

# This loop assigns random initial phases on the interval [0,2*pi]
for i in range(0,num_oscillators):
    varname=str('x'+str(i))
    varnamedict.append(varname)
    kuramotoDS.ics[varname]=random.uniform(0,2*np.pi)

kuramotoDS.name='KuramotoOscillator'

Kuramoto=Generator.Vode_ODEsystem(kuramotoDS)
traj=Kuramoto.compute('test')
pts=traj.sample(dt=1)
f,ax=plt.subplots(figsize=(6,6),dpi=100)
ax.set_xlim([-1.1,1.1])
ax.set_ylim([-1.1,1.1])
line,=ax.plot([],[],'o')

Ox=[[]]
Oy=[[]]
time=list(pts['t'])
for i in range(0,len(time)):
    x_summ=[]
    y_summ=[]
    for name in varnamedict:
        x_summ.append(np.cos(pts[name][i]))
        y_summ.append(np.sin(pts[name][i]))
    Oy.append(y_summ)
    Ox.append(x_summ)  
def animate(i):
    line.set_xdata(Ox[i][:])
    line.set_ydata(Oy[i][:])
    return line,
#def init(i):
#    line.set_data([],[])
#    return line,
ani1 = animation.FuncAnimation(f, animate, frames=np.arange(0,len(time)),interval=30)#,blit=True)
theta=np.arange(0,2.1*(np.pi),0.1)
plt.plot(np.cos(theta),np.sin(theta),'k', linestyle='dashed')
#ani1.save('./kuramoto_N15_K009.png',writer='imagemagick')
plt.show()
