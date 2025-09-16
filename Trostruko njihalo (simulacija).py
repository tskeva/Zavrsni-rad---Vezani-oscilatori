import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from matplotlib.animation import FuncAnimation

def deg(kut):
    return kut/180*np.pi

# Runge-Kutta 4 za tri diferencijalne jednadžbe
def Runge_Kutta4(a,b,N,dif_jed1,dif_jed2,dif_jed3,poc_uvjeti):
    h=(b-a)/N
    red=len(poc_uvjeti[0])
    argument=poc_uvjeti
    lista_t=[a]
    lista_y=[argument[0][0]]
    lista_y2=[argument[1][0]]
    lista_y3=[argument[2][0]]
    for j in range(N):
        t=a+j*h
        delta=[0]*red
        t_pomocni=t
        argument_pomocni=argument.copy()
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1],argument[2]]
        for i in range(red):
            delta[i]+=h/6*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h/2)
        t_pomocni=t+h/2
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1],argument[2]]
        for i in range(red):
            delta[i]+=h/3*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h/2)
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1],argument[2]]
        for i in range(red):
            delta[i]+=h/3*lista_der[i]
            argument_pomocni[0].append(argument[0][i]+lista_der[i]*h)
        t_pomocni=t+h
        lista_der=argument_pomocni[0][1:]+[dif_jed1(t_pomocni,argument_pomocni)]
        argument_pomocni=[[],argument[1],argument[2]]
        for i in range(red):
            delta[i]+=h/6*lista_der[i]

        t=a+j*h
        delta2=[0]*red
        t_pomocni=t
        argument_pomocni=argument.copy()
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[],argument[2]]
        for i in range(red):
            delta2[i]+=h/6*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h/2)
        t_pomocni=t+h/2
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[],argument[2]]
        for i in range(red):
            delta2[i]+=h/3*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h/2)
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[],argument[2]]
        for i in range(red):
            delta2[i]+=h/3*lista_der[i]
            argument_pomocni[1].append(argument[1][i]+lista_der[i]*h)
        t_pomocni=t+h
        lista_der=argument_pomocni[1][1:]+[dif_jed2(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],[],argument[2]]
        for i in range(red):
            delta2[i]+=h/6*lista_der[i]

        t=a+j*h
        delta3=[0]*red
        t_pomocni=t
        argument_pomocni=argument.copy()
        lista_der=argument_pomocni[2][1:]+[dif_jed3(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],argument[1],[]]
        for i in range(red):
            delta3[i]+=h/6*lista_der[i]
            argument_pomocni[2].append(argument[2][i]+lista_der[i]*h/2)
        t_pomocni=t+h/2
        lista_der=argument_pomocni[2][1:]+[dif_jed3(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],argument[1],[]]
        for i in range(red):
            delta3[i]+=h/3*lista_der[i]
            argument_pomocni[2].append(argument[2][i]+lista_der[i]*h/2)
        lista_der=argument_pomocni[2][1:]+[dif_jed3(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],argument[1],[]]
        for i in range(red):
            delta3[i]+=h/3*lista_der[i]
            argument_pomocni[2].append(argument[2][i]+lista_der[i]*h)
        t_pomocni=t+h
        lista_der=argument_pomocni[2][1:]+[dif_jed3(t_pomocni,argument_pomocni)]
        argument_pomocni=[argument[0],argument[1],[]]
        for i in range(red):
            delta3[i]+=h/6*lista_der[i]

        for i in range(red):
            argument[0][i]+=delta[i]
            argument[1][i]+=delta2[i]
            argument[2][i]+=delta3[i]

        lista_t.append(t+h)
        lista_y.append(float(argument[0][0]))
        lista_y2.append(float(argument[1][0]))
        lista_y3.append(float(argument[2][0]))
    return lista_t,lista_y,lista_y2,lista_y3

# Parametri
m=1
g=9.81
l=2

Y1=lambda t,y:0*t+((2-(np.cos(y[1][0]-y[2][0]))**2)*(-3*g/l*np.sin(y[0][0])-2*(y[1][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[0][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((-3/2*np.cos(y[0][0]-y[1][0])+1/2*np.cos(y[0][0]+y[1][0]-2*y[2][0]))*(-2*g/l*np.sin(y[1][0])+2*(y[0][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((-np.cos(y[0][0]-y[2][0])+np.cos(y[0][0]-2*y[1][0]+y[2][0]))*(-g/l*np.sin(y[2][0])+(y[0][1])**2*np.sin(y[0][0]-y[2][0])+(y[1][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)
Y2=lambda t,y:0*t+((-3/2*np.cos(y[0][0]-y[1][0])+1/2*np.cos(y[0][0]+y[1][0]-2*y[2][0]))*(-3*g/l*np.sin(y[0][0])-2*(y[1][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[0][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((3-(np.cos(y[0][0]-y[2][0]))**2)*(-2*g/l*np.sin(y[1][0])+2*(y[0][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((-2*np.cos(y[1][0]-y[2][0])+np.cos(-2*y[0][0]+y[1][0]+y[2][0]))*(-g/l*np.sin(y[2][0])+(y[0][1])**2*np.sin(y[0][0]-y[2][0])+(y[1][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)
Y3=lambda t,y:0*t+((-np.cos(y[0][0]-y[2][0])+np.cos(y[0][0]-2*y[1][0]+y[2][0]))*(-3*g/l*np.sin(y[0][0])-2*(y[1][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[0][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((-2*np.cos(y[1][0]-y[2][0])+np.cos(-2*y[0][0]+y[1][0]+y[2][0]))*(-2*g/l*np.sin(y[1][0])+2*(y[0][1])**2*np.sin(y[0][0]-y[1][0])-(y[2][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)+((4-2*np.cos(2*y[0][0]-2*y[1][0]))*(-g/l*np.sin(y[2][0])+(y[0][1])**2*np.sin(y[0][0]-y[2][0])+(y[1][1])**2*np.sin(y[1][0]-y[2][0])))/(-np.cos(2*y[0][0]-2*y[1][0])-1/2*np.cos(2*y[1][0]-2*y[2][0])+5/2)

# Primjer početnih uvjeta i integracije
t1,y1,y2,y3=Runge_Kutta4(0,10,10000,Y1,Y2,Y3,[[deg(30),-2],[deg(15),1],[deg(-15),-1]])


plt.plot(t1,y1,'r',label='Prvo njihalo')
plt.plot(t1,y2,'b',label='Drugo njihalo')
plt.plot(t1,y3,'g',label='Treće njihalo')
plt.xlabel('t[s]')
plt.ylabel('θ[rad]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=30°, ω\u2081=-2, θ\u2082=15°, ω\u2082=1, θ\u2083=-15°, ω\u2083=-1')
plt.title('θ-t graf')
plt.legend()
plt.show()

X1=[]
Y1=[]

for i in range(len(t1)):
    x=l*np.sin(y1[i])
    y=-l*np.cos(y1[i])
    X1.append(x)
    Y1.append(y)

X2=[]
Y2=[]

for i in range(len(t1)):
    x=l*np.sin(y1[i])+l*np.sin(y2[i])
    y=-l*np.cos(y1[i])-l*np.cos(y2[i])
    X2.append(x)
    Y2.append(y)

X3=[]
Y3=[]

for i in range(len(t1)):
    x=l*np.sin(y1[i])+l*np.sin(y2[i])+l*np.sin(y3[i])
    y=-l*np.cos(y1[i])-l*np.cos(y2[i])-l*np.cos(y3[i])
    X3.append(x)
    Y3.append(y)


plt.plot(t1,Y1,'r',label='Prvo njihalo')
plt.plot(t1,Y2,'b',label='Drugo njihalo')
plt.plot(t1,Y3,'g',label='Treće njihalo')
plt.xlabel('t[s]')
plt.ylabel('y[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=30°, ω\u2081=-2, θ\u2082=15°, ω\u2082=1, θ\u2083=-15°, ω\u2083=-1')
plt.title('y-t graf')
plt.legend()
plt.show()
plt.plot(t1,X1,'r',label='Prvo njihalo')
plt.plot(t1,X2,'b',label='Drugo njihalo')
plt.plot(t1,X3,'g',label='Treće njihalo')
plt.xlabel('t[s]')
plt.ylabel('x[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=30°, ω\u2081=-2, θ\u2082=15°, ω\u2082=1, θ\u2083=-15°, ω\u2083=-1')
plt.title('x-t graf')
plt.legend()
plt.show()
plt.plot(X1,Y1,'r',label='Prvo njihalo')
plt.plot(X2,Y2,'b',label='Drugo njihalo')
plt.plot(X3,Y3,'g',label='Treće njihalo')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=30°, ω\u2081=-2, θ\u2082=15°, ω\u2082=1, θ\u2083=-15°, ω\u2083=-1')
plt.title('y-x graf')
plt.legend()
plt.show()

# Koordinate
X1=[l*np.sin(y1[i]) for i in range(len(t1))]
Yc1=[-l*np.cos(y1[i]) for i in range(len(t1))]

X2=[l*np.sin(y1[i])+l*np.sin(y2[i]) for i in range(len(t1))]
Yc2=[-l*np.cos(y1[i])-l*np.cos(y2[i]) for i in range(len(t1))]

X3=[l*np.sin(y1[i])+l*np.sin(y2[i])+l*np.sin(y3[i]) for i in range(len(t1))]
Yc3=[-l*np.cos(y1[i])-l*np.cos(y2[i])-l*np.cos(y3[i]) for i in range(len(t1))]

# Animacija trostrukog njihala
trajanje_animacije = t1[-1] - t1[0]
frame_rate = 60
broj_frameova = int(trajanje_animacije * frame_rate)

indeksi = np.linspace(0, len(t1) - 1, broj_frameova, dtype=int)
t_anim = [t1[i] for i in indeksi]
X1_anim = [X1[i] for i in indeksi]
Y1_anim = [Yc1[i] for i in indeksi]
X2_anim = [X2[i] for i in indeksi]
Y2_anim = [Yc2[i] for i in indeksi]
X3_anim = [X3[i] for i in indeksi]
Y3_anim = [Yc3[i] for i in indeksi]

fig, ax = plt.subplots(1, 1, figsize=(6,6))
ax.set_facecolor("w")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
ax.set_title("Simulacija trostrukog njihala")

ln_path, = ax.plot([], [], 'k-', lw=1.5, label='Time: 0')
ln1, = ax.plot([], [], 'ro-', lw=2, markersize=9, label='Prvo njihalo')
ln2, = ax.plot([], [], 'bo-', lw=2, markersize=9, label='Drugo njihalo')
ln3, = ax.plot([], [], 'go-', lw=2, markersize=9, label='Treće njihalo')
L = ax.legend(loc=1)

ax.set_ylim(-3*l, l)
ax.set_xlim(-3*l, 3*l)

def animate(i):
    ln_path.set_data([0, X1_anim[i], X2_anim[i], X3_anim[i]],
                     [0, Y1_anim[i], Y2_anim[i], Y3_anim[i]])
    ln1.set_data([X1_anim[i]], [Y1_anim[i]])
    ln2.set_data([X2_anim[i]], [Y2_anim[i]])
    ln3.set_data([X3_anim[i]], [Y3_anim[i]])
    lab = f'Time: {t_anim[i]:.2f} s'
    L.get_texts()[0].set_text(lab)
    return ln_path, ln1, ln2, ln3, L

ani = animation.FuncAnimation(fig, animate, frames=broj_frameova,
                              interval=1000/frame_rate, blit=False, repeat=False)

plt.show()
