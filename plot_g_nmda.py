import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('signal.dat');
t = data[:,0]
V = data[:,1]
g_nmda_fast = data[:,2]
g_nmda_slow = data[:,3]
g_nmda = data[:,4]
I_syn_nmda = data[:,5]
I_Ca_nmda = data[:,6]
I_Na_kdr = data[:,7]

plt.figure(1)
plt.plot(t, V,label="Spine compartment potential",linewidth=2.0)
plt.xlabel("Time [ms]")
plt.ylabel("Spine compartment potential [mV]")
plt.legend()
#plt.figure(2)
#plt.plot(t, g_nmda_fast)
#plt.figure(3)
#plt.plot(t, g_nmda_slow)
#plt.figure(4)
#plt.plot(t, g_nmda)
#plt.figure(5)
#plt.plot(t, I_syn_nmda)
#plt.figure(6)
#plt.plot(t, I_Ca_nmda)
#plt.figure(7)
#plt.plot(t, I_Na_kdr)
plt.show();
