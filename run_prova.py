import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import nest
from scipy.optimize import curve_fit
import numpy as np
import os
from scipy.signal import find_peaks

from pynestml.codegeneration.nest_code_generator_utils import NESTCodeGeneratorUtils

# generate and build code
module_name, neuron_model_name = \
    NESTCodeGeneratorUtils.generate_code_for("hh_cond_exp_Ca_bro.nestml")



import nest
nest.set_verbosity("M_ALL")
nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.01, "print_time": True} )
nest.Install(module_name)
font_lb=14
neuron = nest.Create(neuron_model_name)
voltmeter = nest.Create("voltmeter")
voltmeter.set({"record_from": ["V_spine", "I_Ca_spine","m_ca_sp","g_nmda","h_ca_sp","Ca_spine" ,"Ca_spine1" ,"s_na_sp", "I_Ca_spine_NMDA","I_Ca_spine_NMDA1","I_Na_spine","I_mAHP_spine","h_na_sp","I_Kdr_spine","m_kdr_sp","m_na_soma","h_na_soma","s_na_soma","I_Na_soma","m_mAHP_sp","m_a_sp","h_a_sp","I_A_spine","I_leak","g_nmda_fast","g_nmda_slow","I_syn_spine_NMDA","h_nmda","I_syn_exc","m_na_sp","h_na_sp_prova","m_na_sp_prova","m_kdr_prova"],"interval": 0.01})
#neuron.SetStatus(neuron,{"I_e":10})
nest.Connect(voltmeter, neuron)


#sg = nest.Create('spike_generator', params={"spike_times": [0.01]})
#nest.Connect(sg, neuron, syn_spec={"weight": 1.0, "delay": 0.1})
nest.Simulate(300.)


data = np.loadtxt('signal.dat')
#data_1 = np.loadtxt('currents_signals.dat')
t_cpp = data[:,0]
V_cpp = data[:,1]
g_nmda_fast_cpp=data[:,2]
g_nmda_slow_cpp=data[:,3]
g_nmda_cpp=data[:,4]
i_syn_cpp=data[:,5]
i_ca_cpp=data[:,6]
i_na_cpp=data[:,7]
i_k_cpp=data[:,8]
i_l_cpp=data[:,9]
m_cpp=data[:,10]
h_cpp=data[:,11]
n_cpp=data[:,12]
s_cpp=data[:,13]
#i_l = data[:,]
#i_na = data[:,3]
#i_kdr = data[:,4]

t=voltmeter.get("events")["times"]
V_spine_nest=voltmeter.get("events")["V_spine"]
I_L_nest=voltmeter.get("events")["I_leak"]
I_Na_nest=voltmeter.get("events")["I_Na_spine"]
ca_concentration= voltmeter.get("events")["Ca_spine"]
i_kdr_nest=voltmeter.get("events")["I_Kdr_spine"]
i_nmda_nest=voltmeter.get("events")["I_syn_spine_NMDA"]
i_ca_nmda_nest=voltmeter.get("events")["I_Ca_spine_NMDA"]
h_nmda_nest=voltmeter.get("events")["h_nmda"]
g_nmda_nest=voltmeter.get("events")["g_nmda"]
i_syn_exc=voltmeter.get("events")["I_syn_exc"]
m_nest=voltmeter.get("events")["m_na_sp_prova"]
h_nest=voltmeter.get("events")["h_na_sp_prova"]
#n_nest=voltmeter.get("events")["m_kdr_sp"]
n_nest=voltmeter.get("events")["m_kdr_prova"]
s_nest=voltmeter.get("events")["s_na_sp"]

g_nmda_fast_nest = voltmeter.get("events")["g_nmda_fast"]
g_nmda_slow_nest = voltmeter.get("events")["g_nmda_slow"]

fig2, ax2 = plt.subplots(7) #nrows=3,ncols=1)
ax2[0].plot(t_cpp,V_cpp, label="C++",linewidth=5.0)
ax2[0].plot(t,V_spine_nest, label="NESTML",linewidth=5.0,linestyle="--")
ax2[0].set_ylabel(r"Spine potential",fontsize=font_lb)

#ax1[0].plot(t, g_nmda_nest, label="NESTML",linewidth=2.0,linestyle="-")
#ax1[0].set_ylabel(r"g_nmda",fontsize=font_lb)
ax2[0].set_xlabel(r"Time [ms]",fontsize=font_lb)
#ax1[0][0].set_xlim([0,10])
ax2[0].legend()



ax2[1].plot(t, i_nmda_nest/1e6, label="NESTML",linewidth=5.0,linestyle="-")
ax2[1].plot(t_cpp, i_syn_cpp, label="C++",linewidth=5.0,linestyle="--")
ax2[1].set_ylabel(r"i_nmda",fontsize=font_lb)
ax2[1].set_xlabel(r"Time [ms]",fontsize=font_lb)
#a1[0][0].set_xlim([0,10])
ax2[1].legend()

#ax1[2].plot(t_cpp,i_na, label="C++",linewidth=5.0)
#ax1[2].plot(t,I_Na_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
#ax1[2].set_ylabel(r"Sodium current",fontsize=font_lb)
#ax1[2].set_xlabel(r"Time [ms]",fontsize=font_lb)

ax2[2].plot(t, i_ca_nmda_nest/1e6, label="NESTML",linewidth=5.0,linestyle="-")
ax2[2].plot(t_cpp, i_ca_cpp, label="C++",linewidth=5.0,linestyle="--")
ax2[2].set_ylabel(r"i_ca_nmda",fontsize=font_lb)
#ax1[0][0].set_xlim([0,10])
ax2[2].legend()

#ax1[3].plot(t_cpp,i_kdr, label="C++",linewidth=5.0)
#ax1[3].plot(t,i_kdr_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
#ax1[3].set_ylabel(r"Potassium current",fontsize=font_lb)
#ax1[3].set_xlabel(r"Time [ms]",fontsize=font_lb)
#ax1[0][0].set_xlim([0,10])
#ax1[3].legend()

#ax1[4].plot(t,i_nmda_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
#ax1[4].set_ylabel(r"I syn nmda",fontsize=font_lb)
#ax1[4].set_xlabel(r"Time [ms]",fontsize=font_lb)

ax2[3].plot(t_cpp,-i_na_cpp, label="C++",linewidth=5.0)
ax2[3].plot(t,I_Na_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
ax2[3].set_ylabel(r"ina",fontsize=font_lb)
ax2[3].legend()

ax2[4].plot(t_cpp,-i_k_cpp, label="C++",linewidth=5.0)
ax2[4].plot(t,i_kdr_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
ax2[4].set_ylabel(r"ik",fontsize=font_lb)
ax2[4].legend()

ax2[5].plot(t_cpp,-i_l_cpp, label="C++",linewidth=5.0)
ax2[5].plot(t,I_L_nest/1e6, label="NESTML",linewidth=5.0,linestyle="--")
ax2[5].set_ylabel(r"il",fontsize=font_lb)
ax2[5].legend()

ax2[6].plot(t,ca_concentration)
plt.xticks([])
indici_picchi, _ = find_peaks(V_spine_nest, height=23.0)
plt.xticks(t[indici_picchi], range(1, len(indici_picchi) + 1)) 
ax2[6].set_ylabel(r"Calcium concentration",fontsize=font_lb)
ax2[6].legend()

fig1, ax1 = plt.subplots(3)
ax1[0].plot(t_cpp,m_cpp, label="C++",linewidth=5.0)
ax1[0].plot(t,m_nest, label="NESTML",linewidth=5.0,linestyle="--")
ax1[0].set_ylabel(r"m sodium",fontsize=font_lb)
ax1[0].legend()

ax1[1].plot(t_cpp,h_cpp, label="C++",linewidth=5.0)
ax1[1].plot(t,h_nest, label="NESTML",linewidth=5.0,linestyle="--")
ax1[1].set_ylabel(r"h sodium",fontsize=font_lb)
ax1[1].legend()


ax1[2].plot(t_cpp,n_cpp, label="C++",linewidth=5.0)
ax1[2].plot(t,n_nest, label="NESTML",linewidth=5.0,linestyle="--")
ax1[2].set_ylabel(r"n pot",fontsize=font_lb)
ax1[2].legend()


#ax1[3].plot(t_cpp,s_cpp, label="C++",linewidth=5.0)
#ax1[3].plot(t,s_nest, label="NESTML",linewidth=5.0,linestyle="--")
#ax1[3].set_ylabel(r"s",fontsize=font_lb)
#ax1[3].legend()

#plt.plot(t,i_syn_exc)



plt.figure(3)
plt.plot(t,ca_concentration)
plt.xticks([])
indici_picchi, _ = find_peaks(V_spine_nest, height=23.0)
plt.xticks(t[indici_picchi], range(1, len(indici_picchi) + 1)) 
plt.ylabel(r"[$Ca^{2+}$]",fontsize=font_lb)
plt.xlabel(r"Action potential #",fontsize=font_lb)

plt.legend()
plt.plot()

plt.show()

