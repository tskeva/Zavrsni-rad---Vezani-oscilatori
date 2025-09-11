import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation


def deg(kut):
    return kut/180*np.pi

#Rješavanje sustava N diferencijalnih jednadžbi N-tog reda

def Runge_Kutta4(a,b,N,dif_jed1,dif_jed2,poc_uvjeti):
    '''dif_jedn: yN=f(t,[y0,...,yN-1])\npoc_uvjeti: [y0(t0),...,yN-1(t0)]'''
    h=(b-a)/N
    red=len(poc_uvjeti)
    argument=poc_uvjeti
    lista_t=[a]
    lista_y=[argument[0][0]]
    lista_y2=[argument[1][0]]
    for j in range(N):
        t=a+j*h
        delta=[0]*red
        t_pomocni=t
        argument_pomocni=argument.copy()
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1]]
        for i in range(red):
            delta[i]+=h/6*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h/2)
        t_pomocni=t+h/2
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1]]
        for i in range(red):
            delta[i]+=h/3*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h/2)
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1]]
        for i in range(red):
            delta[i]+=h/3*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h)
        t_pomocni=t+h
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1]]
        for i in range(red):
            delta[i]+=h/6*lista_der[i]

        t=a+j*h
        delta2=[0]*red
        t_pomocni=t
        argument_pomocni=argument.copy()
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[]]
        for i in range(red):
            delta2[i]+=h/6*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h/2)
        t_pomocni=t+h/2
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[]]
        for i in range(red):
            delta2[i]+=h/3*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h/2)
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[]]
        for i in range(red):
            delta2[i]+=h/3*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h)
        t_pomocni=t+h
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[]]
        for i in range(red):
            delta2[i]+=h/6*lista_der[i]
            argument[0][i]+=delta[i]
            argument[1][i]+=delta2[i]
        lista_t.append(t+h)
        lista_y.append(float(argument[0][0]))
        lista_y2.append(float(argument[1][0]))
    return lista_t,lista_y,lista_y2


m1=0.2 #kg
m2=0.1 #kg
l1=1.25 #m
l2=0.75 #m

Y1=lambda t,y:0*t+(-9.81*(2*m1+m2)*np.sin(y[0][0])-m2*9.81*np.sin(y[0][0]-2*y[1][0])-2*np.sin(y[0][0]-y[1][0])*m2*((y[1][1])**2*l2+(y[0][1])**2*l1*np.cos(y[0][0]-y[1][0])))/(l1*(2*m1+m2-m2*np.cos(2*y[0][0]-2*y[0][1])))
Y2=lambda t,y:0*t+(2*np.sin(y[0][0]-y[1][0])*((y[0][1])**2*l1*(m1+m2)+9.81*(m1+m2)*np.cos(y[0][0])+(y[1][1])**2*l2*m2*np.cos(y[0][0]-y[1][0])))/(l2*(2*m1+m2-m2*np.cos(2*y[0][0]-2*y[1][0])))

t1,y1,y2=Runge_Kutta4(0,15,10000,Y1,Y2,[[deg(15),-2],[deg(30),4]]) #Potrebno unijeti: t1,t2, broj koraka N, prva dif. jed., druga dif. jed., poc. uvjeti za oba njihala

#Grafovi

plt.plot(t1,y1,'r',label='Prvo njihalo')
plt.plot(t1,y2,'b',label='Drugo njihalo')
plt.xlabel('t[s]')
plt.ylabel('θ[rad]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=15°, ω\u2081=-2, θ\u2082=30°, ω\u2082=4')
plt.title('θ-t graf')
plt.legend()
plt.show()

X1=[]
Y1=[]

for i in range(len(t1)):
    x=l1*np.sin(y1[i])
    y=-l1*np.cos(y1[i])
    X1.append(x)
    Y1.append(y)

X2=[]
Y2=[]

for i in range(len(t1)):
    x=l1*np.sin(y1[i])+l2*np.sin(y2[i])
    y=-l1*np.cos(y1[i])-l2*np.cos(y2[i])
    X2.append(x)
    Y2.append(y)


plt.plot(t1,Y1,'r',label='Prvo njihalo')
plt.plot(t1,Y2,'b',label='Drugo njihalo')
plt.xlabel('t[s]')
plt.ylabel('y[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=15°, ω\u2081=-2, θ\u2082=30°, ω\u2082=4')
plt.title('y-t graf')
plt.legend()
plt.show()
plt.plot(t1,X1,'r',label='Prvo njihalo')
plt.plot(t1,X2,'b',label='Drugo njihalo')
plt.xlabel('t[s]')
plt.ylabel('x[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=15°, ω\u2081=-2, θ\u2082=30°, ω\u2082=4')
plt.title('x-t graf')
plt.legend()
plt.show()
plt.plot(X1,Y1,'r',label='Prvo njihalo')
plt.plot(X2,Y2,'b',label='Drugo njihalo')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=15°, ω\u2081=-2, θ\u2082=30°, ω\u2082=4')
plt.title('y-x graf')
plt.legend()
plt.show()

#Animacija dvostrukog njihala

trajanje_animacije = t1[-1] - t1[0]
frame_rate = 60
broj_frameova = int(trajanje_animacije * frame_rate)

indeksi = np.linspace(0, len(t1) - 1, broj_frameova, dtype=int)
t_anim = [t1[i] for i in indeksi]
X1_anim = [X1[i] for i in indeksi]
Y1_anim = [Y1[i] for i in indeksi]
X2_anim = [X2[i] for i in indeksi]
Y2_anim = [Y2[i] for i in indeksi]

fig, ax = plt.subplots(1, 1, figsize=(6,6))
ax.set_facecolor("w")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
ax.set_title("Simulacija vezanih njihala")

ln4, = ax.plot([], [], 'k-', lw=1.5, label='Time: 0')
ln1, = ax.plot([], [], 'k')
ln2, = ax.plot([], [], 'ro-', lw=2, markersize=9, label='Prvo njihalo')
ln3, = ax.plot([], [], 'bo-', lw=2, markersize=9, label='Drugo njihalo')
L = ax.legend(loc=1)

ax.set_ylim(-4, 1)
ax.set_xlim(-2.5, 2.5)

def animate(i):
    ln1.set_data([0, X1_anim[i], X2_anim[i]], [0, Y1_anim[i], Y2_anim[i]])
    ln2.set_data([X1_anim[i]], [Y1_anim[i]])
    ln3.set_data([X2_anim[i]], [Y2_anim[i]])
    ln4.set_data([0, X1_anim[i], X2_anim[i]], [0, Y1_anim[i], Y2_anim[i]])
    lab = f'Time: {t_anim[i]:.2f} s'
    L.get_texts()[0].set_text(lab)
    return ln1, ln2, ln3, ln4, L

ani = animation.FuncAnimation(fig, animate,frames=broj_frameova,interval=1000 / frame_rate,blit=False,repeat=False)

plt.show()
