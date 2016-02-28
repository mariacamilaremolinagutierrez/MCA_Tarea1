import os
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True, fontsize=22)

#x, rho, u, e, P
data = np.loadtxt('final_step.dat')
data_lf = np.loadtxt('final_step_laxfriedrichs.dat')

os.system('mkdir Graficas/')

x = data[:,0]
rho = data[:,1]
u = data[:,2]
e = data[:,3]
P = data[:,4]

x_lf = data_lf[:,0]
rho_lf = data_lf[:,1]
u_lf = data_lf[:,2]
e_lf = data_lf[:,3]
P_lf = data_lf[:,4]

fig = plt.figure(figsize=(16,8))
plt.plot(x,rho,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_lf,rho_lf,label="$\mathrm{Lax-Friedrichs\'s\;Scheme}$")
plt.ylim(-0.1,1.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Density}$')
plt.legend(fontsize=20)
plt.savefig('Graficas/density.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,u,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_lf,u_lf,label="$\mathrm{Lax-Friedrichs\'s\;Scheme}$")
plt.ylim(-0.1,0.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Velocity}$')
plt.legend(fontsize=20)
plt.savefig('Graficas/velocity.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,e,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_lf,e_lf,label="$\mathrm{Lax-Friedrichs\'s\;Scheme}$")
plt.ylim(-0,2.7)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Energy\;per\;unit\;lenghth}$')
plt.legend(fontsize=20)
plt.savefig('Graficas/energyperunitlength.png')

fig = plt.figure(figsize=(16,8))
plt.plot(x,P,label="$\mathrm{Godunov\'s\;Scheme}$")
plt.plot(x_lf,P_lf,label="$\mathrm{Lax-Friedrichs\'s\;Scheme}$")
plt.ylim(-0.1,1.1)
plt.xlabel('$\mathrm{x}$')
plt.ylabel('$\mathrm{Pressure}$')
plt.legend(fontsize=20)
plt.savefig('Graficas/pressure.png')

print np.average(rho_lf)
