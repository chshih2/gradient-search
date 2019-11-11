import math
import numpy as np
from scipy.integrate import ode
from numpy.linalg import inv
import os
import glob
import shutil
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
# theta1_desired=2.5
# theta2_desired=-2.7

# targetx=0.13
# targety=0.05
theta1_0=3.0
theta2_0=-2.0

# coef_v=0.002

F1_limit=0.08
F2_limit=0.12

density = 1000.0 #rod density
L0 = 0.1 #rod length
r0 = 0.01 #rod radius
totalMass = density*math.pi*r0*r0*L0 #total mass for one rod
u1=0.01 #scalar torque at joint 1
u2=0.00 #scalar torque at joint 2

# %% Robot Specifications
M1=totalMass; #link 1 mass
M2=totalMass; #link 2 mass
M=totalMass;
L1=L0; #link 1 length
L2=L0; #link 1 length
# L=L0;
def LnU(limit,c,b,m,s_dot):
    Fmin=-limit
    Fmax=limit
    dummy_a=(Fmin-c*s_dot**2-b*s_dot)/m
    dummy_b=(Fmax-c*s_dot**2-b*s_dot)/m
    if m>0:
        L=dummy_a
        U=dummy_b
    elif m<0:
        L=dummy_b
        U=dummy_a
    else:
        print("m=0")
    return L,U
def theta_s(s,s_dot):
    theta1_s=theta1_0+s*(theta1_desired-theta1_0)
    theta2_s=theta2_0+s*(theta2_desired-theta2_0)
    m1_s=M*L0*L0*((3+2*np.cos(theta2_s))*dtheta1_ds+(3+np.cos(theta2_s))*dtheta2_ds)
    b1_s=coef_v*dtheta1_ds
    c1_s=-M*L0*L0*np.sin(theta2_s)*(2*dtheta1_ds*dtheta2_ds+dtheta2_ds*dtheta2_ds)
    m2_s=M*L0*L0*((1+np.cos(theta2_s))*dtheta1_ds+dtheta2_ds)
    b2_s=coef_v*dtheta2_ds
    c2_s=-M*L0*L0*np.sin(theta2_s)*(dtheta1_ds*dtheta1_ds)
    L1,U1=LnU(F1_limit,c1_s,b1_s,m1_s,s_dot)
    L2,U2=LnU(F2_limit,c2_s,b2_s,m2_s,s_dot)
    L=np.maximum(L1,L2)
    U=np.minimum(U1,U2)
    # print(str(L)+"  "+str(U))
    return L,U,theta1_s,theta2_s
# def robot(t, q,acc): #lumped system dynamics https://ocw.mit.edu/courses/mechanical-engineering/2-12-introduction-to-robotics-fall-2005/lecture-notes/chapter7.pdf
#     s=q[0]
#     s_dot=acc*dt
#     return [s_dot,acc]



coef_v=0.002 #viscous joint friction coefficient
coef_c=0.00 #coulomb joint friction coefficient

def inversekinematics(flag):
    x=targetx
    y=targety
    alpha=math.acos((x**2+y**2+L1**2-L2**2)/(2*L1*np.sqrt(x**2+y**2)))

    beta=math.acos((L1**2+L2**2-x**2-y**2)/(2*L1*L2))
    gamma=math.atan2(x,y)
    if flag==0:
        theta1_desired=gamma+alpha
        theta2_desired=beta-math.pi
    else:
        theta1_desired=gamma-alpha
        theta2_desired=math.pi-beta
    return theta1_desired,theta2_desired
record_flag=[]
# for j in range(100):
for j in [1]:
    # length = np.sqrt(np.random.uniform(0, 1))*0.2
    # angle = np.pi * np.random.uniform(0, 2)
    # targetx=length * np.cos(angle)
    # targety=length * np.sin(angle)
    targetx=0.09104983281833028
    targety=-0.01561079120308262
    flag=0
    while flag==0 or flag==1:

        theta1_desired,theta2_desired=inversekinematics(flag)
        dtheta1_ds=(theta1_desired-theta1_0)
        dtheta2_ds=(theta2_desired-theta2_0)
        time=0.001
        # ================backward================
        dt=-time #time step
        backward=[]
        s=1
        s_dot=0
        s_dot_old=0
        acc=0
        t=0

        while True:
            s = s+s_dot*dt
            s_dot = s_dot+acc*dt

            backward.append((s,s_dot))

            L,U,_,_=theta_s(s,s_dot)
            if (U>L and s>0):
                acc=L
                if s_dot<s_dot_old:
                    break
            else:
                break
            s_dot_old=s_dot
        x,y=np.array(backward).T
        plt.plot(x,y,'ro-')
        Lower=L
        # ================forward================
        dt=time #time step
        forward=[]
        s=0
        s_dot=0
        acc=0
        t=0
        goal=0
        count=0
        print("weird "+str(j)+" target: "+str(targetx)+"    "+str(targety))

        while True:
            s = s+s_dot*dt
            s_dot = s_dot+acc*dt
            print(str(s)+"  "+str(s_dot))
            forward.append((s,s_dot))

            L,U,_,_=theta_s(s,s_dot)
            if U>L:
                acc=U
                count+=1
                if np.min(cdist([(s,s_dot)],backward))<0.03 or s>1:
                    goal=1
                    break
            else:
                flag+=1
                break
            # print("dis  "+str())
        x,y=np.array(forward).T
        plt.plot(x,y,'bo-')
        Upper=U
        # plt.show()
        plt.savefig("maxminLnU"+str(j)+".png")
        plt.close()

        if goal==1:
            break
    record_flag.append(goal)
print(record_flag)
print(len(np.where(record_flag==1)))
x1_desired=L1*np.sin(theta1_desired) #calculate position
y1_desired=L1*np.cos(theta1_desired)
x2_desired=x1_desired+L2*np.sin(theta1_desired+theta2_desired)
y2_desired=y1_desired+L2*np.cos(theta1_desired+theta2_desired)

# ================backward================
dt=-time #time step
backward=[]
s=1
s_dot=0
acc=Lower
t=0
while True:
    s = s+s_dot*dt
    s_dot = s_dot+acc*dt

    backward.append((s,s_dot))
    if s<0 or s>1:
        print("oh")
        break
x,y=np.array(backward).T
plt.plot(x,y,'ro-')
# ================forward================
dt=time #time step
forward=[]
s=0
s_dot=0
acc=Upper
t=0
while True:
    s = s+s_dot*dt
    s_dot = s_dot+acc*dt

    forward.append((s,s_dot))

    print(s)
    if s<0 or s>1 or np.min(cdist([(s,s_dot)],backward))<0.03:
        break
x,y=np.array(forward).T
plt.plot(x,y,'bo-')
Upper=U
plt.savefig("phase.png")
plt.close()
switch=(x[-1],y[-1])
print(switch)
#
theta1_record=[]
theta2_record=[]
phase_record=[]
dt=time #time step
s=0
s_dot=0
acc=Upper
T=0
while True:
    s = s+s_dot*dt
    s_dot = s_dot+acc*dt
    phase_record.append((s,s_dot))
    _,_,theta1,theta2=theta_s(s,s_dot)
    theta1_record.append(theta1)
    theta2_record.append(theta2)

    x1=L1*np.sin(theta1) #calculate position
    y1=L1*np.cos(theta1)
    x2=x1+L2*np.sin(theta1+theta2)
    y2=y1+L2*np.cos(theta1+theta2)


    if switch[0]<s:
        acc=Lower

    if np.min(cdist([(x2,y2)],[(x2_desired,y2_desired)]))<0.01:
        break
    T+=dt
print(T)
x,y=np.array(phase_record).T
plt.plot(x,y,'ro-')
plt.savefig("phase_portrait.png")
plt.close()
# print(len(theta1_record))
#
totalframe=list(range(0,len(theta1_record),10))
totalframe.append(len(theta1_record)-1)
totalframe.append(len(theta1_record)-1)
print(totalframe)

for i in range(0, len(totalframe)):
    theta1=theta1_record[totalframe[i]]
    theta2=theta2_record[totalframe[i]]
    x1=L1*np.sin(theta1) #calculate position
    y1=L1*np.cos(theta1)
    x2=x1+L2*np.sin(theta1+theta2)
    y2=y1+L2*np.cos(theta1+theta2)
    plt.plot([0,x1,x2], [0,y1,y2], 'ro-')#,alpha=0.2*(i/(len(theta1_record)+1)))
    plt.plot(targetx,targety,'ko')
    plt.plot([0,x1_desired,x2_desired], [0,y1_desired,y2_desired], 'bo-',alpha=0.2)
    plt.ylim(-0.25,0.25)
    plt.xlim(-0.25,0.25)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("{:04d}".format(i)+'.png')
    plt.close()
print(len(totalframe)/T)
os.system('ffmpeg -threads 8 -r '+str(len(theta1_record)/T/10)+' -i ./%04d.png -b:v 90M -vsync 2 -vcodec mpeg4 ./0.mp4')
for fl in glob.glob("./0*.png"):
    os.remove(fl)
