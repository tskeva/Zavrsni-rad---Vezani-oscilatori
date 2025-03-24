import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import animation
from matplotlib.animation import PillowWriter
from matplotlib.animation import FuncAnimation


def deg(kut):
    return kut/180*np.pi

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
        #print(type(t+h),type(float(argument[0][0])))


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


m1=0.2
m2=0.1
l1=1.25
l2=0.75

Y1=lambda t,y:0*t+(-9.81*(2*m1+m2)*np.sin(y[0][0])-m2*9.81*np.sin(y[0][0]-2*y[1][0])-2*np.sin(y[0][0]-y[1][0])*m2*((y[1][1])**2*l2+(y[0][1])**2*l1*np.cos(y[0][0]-y[1][0])))/(l1*(2*m1+m2-m2*np.cos(2*y[0][0]-2*y[0][1])))
Y2=lambda t,y:0*t+(2*np.sin(y[0][0]-y[1][0])*((y[0][1])**2*l1*(m1+m2)+9.81*(m1+m2)*np.cos(y[0][0])+(y[1][1])**2*l2*m2*np.cos(y[0][0]-y[1][0])))/(l2*(2*m1+m2-m2*np.cos(2*y[0][0]-2*y[1][0])))

t1,y1,y2=Runge_Kutta4(0,10,10000,Y1,Y2,[[deg(45),2*np.pi/180],[deg(30),4*np.pi/180]]) #0,1000,10000,15,2,30,4
#0,10,1000,15,2,30,4
#0,10,1000,45,2,30,4

#Sustav s u i v varijablama

#print(t1[0],y1[0])
#print(type(t1[0]),type(y1[0]))
plt.plot(t1,y1,'r',label='Prvo njihalo')
plt.plot(t1,y2,'b',label='Drugo njihalo')
plt.xlabel('t[s]')
plt.ylabel('θ[rad]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=45°, ω\u2081=2, θ\u2082=30°, ω\u2082=4')
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
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=45°, ω\u2081=2, θ\u2082=30°, ω\u2082=4')
plt.title('y-t graf')
plt.legend()
plt.show()
plt.plot(t1,X1,'r',label='Prvo njihalo')
plt.plot(t1,X2,'b',label='Drugo njihalo')
plt.xlabel('t[s]')
plt.ylabel('x[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=45°, ω\u2081=2, θ\u2082=30°, ω\u2082=4')
plt.title('x-t graf')
plt.legend()
plt.show()
plt.plot(X1,Y1,'r',label='Prvo njihalo')
plt.plot(X2,Y2,'b',label='Drugo njihalo')
plt.xlabel('x[m]')
plt.ylabel('y[m]')
plt.scatter([], [], color="w", alpha=0, label='θ\u2081=45°, ω\u2081=2, θ\u2082=30°, ω\u2082=4')
plt.title('y-x graf')
plt.legend()
plt.show()


def animate(i):
    ln1.set_data([0,X1[i],X2[i]],[0,Y1[i],Y2[i]])
    ln2.set_data([X1[i]],[Y1[i]])
    ln3.set_data([X2[i]],[Y2[i]])
    ln4.set_data([0,X1[i],X2[i]],[0,Y1[i],Y2[i]])
    lab = 'Time:'+str(round(0.01+0.01*i,1))
    L.get_texts()[0].set_text(lab)

fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.set_facecolor("w")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ln4, =plt.plot([],[], 'k-', lw=1.5,label='Time: 0')
ln1, =plt.plot(0,0,'b')
ln2, =plt.plot([],[], 'ro-', lw=2, markersize=9,label='Prvo njihalo')
ln3, =plt.plot([],[], 'bo-', lw=2, markersize=9,label='Drugo njihalo')
L=plt.legend(loc=1)
ax.set_ylim(-4,1)
ax.set_xlim(-2.5,2.5)
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=0.01, repeat=False)
plt.show()