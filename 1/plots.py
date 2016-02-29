import os
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True, fontsize=22)

#Analytic Solution
#Based on: http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
N = 5000
L = 4.0
dx = L/float(N)
tMax = 1.0
gamma = 1.4
m_2 = (gamma-1)/(gamma + 1)

r_left = 1.0
P_left = 1.0
u_left = 0.0
c_left = np.sqrt(gamma*P_left/r_left)

r_right = 0.125
P_right = 0.1
u_right = 0.0
c_right = np.sqrt(gamma*P_right/r_right)

Pl = P_left
Pr = P_right

while(abs(Pl-Pr) > 0.00001):

    Ps = (Pl-Pr)/2.0 + Pr
    f_Ps = (2*np.sqrt(gamma)/(gamma-1))*(1-Ps**((gamma-1)/(2*gamma))) - (Ps-P_right)*(((1-m_2)**2)/(r_right*(Ps+m_2*P_right)))**0.5

    if(f_Ps < 0):
        Pl = Ps
    else:
        Pr = Ps

P_post = Ps
u_post = (2*np.sqrt(gamma)/(gamma-1))*(1-P_post**((gamma-1)/(2*gamma)))
r_post = r_right*(P_post/P_right+m_2)/(1+m_2*(P_post/P_right))

r_middle = r_left*(P_post/P_left)**(1/gamma)
u_shock = u_post*(r_post/r_right)/(r_post/r_right-1)

x0 = L/2.0
x1 = x0 - c_left*tMax
x2 = x0 + (u_post-(c_left-((gamma-1)/2.0)*u_post))*tMax;
x3 = x0 + u_post*tMax
x4 = x0 + u_shock*tMax

x_a = np.zeros(N)
rho_a = np.zeros(N)
u_a = np.zeros(N)
e_a = np.zeros(N)
P_a = np.zeros(N)

current_x = 0.0

for i in range(0,N):

    if(current_x<x1):
        uu = u_left
        rr = r_left
        PP = P_left

    elif(current_x>x1 and current_x<x2):
        c_sound = m_2*((x0-current_x)/tMax) + (1-m_2)*c_left
        uu = (1-m_2)*(-(x0-current_x)/tMax + c_left)
        rr = r_left*(c_sound/c_left)**(2/(gamma-1))
        PP = P_left*(rr/r_left)**gamma

    elif(current_x>x2 and current_x<x3):
        uu = u_post
        rr = r_middle
        PP = P_post

    elif(current_x>x3 and current_x<x4):
        uu = u_post
        rr = r_post
        PP = P_post

    else:
        uu = u_right
        rr = r_right
        PP = P_right

    ee = PP/((gamma-1)*rr)

    x_a[i] = current_x
    rho_a[i] = rr
    u_a[i] = uu
    P_a[i] = PP
    e_a[i] = ee

    current_x += dx


#Godunov's Solution
#x, rho, u, e, P
data = np.loadtxt('final_step.dat')

os.system('mkdir Graficas/')

x = data[:,0]
rho = data[:,1]
u = data[:,2]
e = data[:,3]
P = data[:,4]

fig = plt.figure(figsize=(16,8))
plt.plot(x,rho,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_a,rho_a,label="$\mathrm{Analytic}$")
plt.ylim(-0.1,1.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Density}$')
plt.legend(fontsize=20, loc='best')
plt.savefig('Graficas/density.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,u,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_a,u_a,label="$\mathrm{Analytic}$")
plt.ylim(-0.1,1.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Velocity}$')
plt.legend(fontsize=20, loc='best')
plt.savefig('Graficas/velocity.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,e,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_a,e_a,label="$\mathrm{Analytic}$")
#plt.ylim(-0,2.7)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Energy\;per\;unit\;lenghth}$')
plt.legend(fontsize=20, loc='best')
plt.savefig('Graficas/energyperunitlength.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,P,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_a,P_a,label="$\mathrm{Analytic}$")
plt.ylim(-0.1,1.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Pressure}$')
plt.legend(fontsize=20, loc='best')
plt.savefig('Graficas/pressure.png')
