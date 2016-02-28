import os
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True, fontsize=22)

#x, rho, u, e, P
data = np.loadtxt('final_step.dat')
#data_a = np.loadtxt('final_step_analytic.dat')
data_a = np.loadtxt('final_step.dat')

os.system('mkdir Graficas/')

x = data[:,0]
rho = data[:,1]
u = data[:,2]
e = data[:,3]
P = data[:,4]

x_a = data_a[:,0]
rho_a = data_a[:,1]
u_a = data_a[:,2]
e_a = data_a[:,3]
P_a = data_a[:,4]

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
plt.ylim(-0,2.7)
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
