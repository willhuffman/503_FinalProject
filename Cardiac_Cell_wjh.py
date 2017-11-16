# Will Huffman
# BME 503 Final Project
# Model based on Dokos 1996 Papers:
# [1]
# "Ion Currents Underlying Sinoatrial Node Pacemaker Activity:
#  A New Single Cell Mathematical Model"
# [2]
# "Vagal Control of Sinoatrial Rhythm] a Mathematical Model"

from brian2 import *
from brian2.units.constants import faraday_constant, gas_constant

defaultclock.dt = 0.01*ms

num_neurons = 1
duration = 10.0*second

eqs_i_CaL = '''
i_CaL = d_L * f_L * f_2L * g_CaL * (v - E_Ca + 75*mV) : amp

d_Linf = 1/(1+exp((v+6.6*mV)/(-6.6*mV))) : 1
tau_d_L = 0.002 * second : second
dd_L/dt   = (d_Linf - d_L)/tau_d_L : 1

f_Linf = 1/(1+exp((v+25*mV)/(6*mV))) : 1
tau_f_L = 0.031*second + 1*second/(1+exp((v+37.6*mV)/(8.1*mV))) : second
df_L/dt = (f_Linf - f_L)/tau_f_L : 1

df_2L/dt = alpha_f_2L * (1 - f_2L) - beta_f_2L*Ca_i*f_2L : 1

g_CaL : siemens
alpha_f_2L : Hz
beta_f_2L : Hz/mM
'''

eqs_i_CaT = '''
i_CaT = d_T * f_T * g_CaT * (v-E_Ca+75*mV) : amp

d_T_inf = 1/(1+exp((v + 23*mV)/(-6.1*mV))) : 1
tau_d_T = 0.0006*second + 0.0054*second/(1+exp(0.03*(v + 100*mV)/mV)) : second
dd_T/dt = (d_T_inf - d_T)/tau_d_T : 1

f_Tinf = 1/(1+exp((v+75*mV)/(6.6*mV))) : 1
tau_f_T = 0.001*second + 0.04*second/(1+exp(0.08*(v+65*mV)/mV)) : second
df_T/dt = (f_Tinf - f_T)/tau_f_T : 1

g_CaT : siemens
'''

eqs_i_Na = '''
i_Na = (m**3) * h * g_Na * (v-E_Na) : amp

dm/dt = alpha_m * (1-m) - beta_m * m : 1

alpha_m = 200*(v + 34.3*mV)/(1-exp(-0.09*(v+34.3*mV)/mV))/(mV*second) : Hz 
beta_m = 8000*exp(-0.15*(v+56.2*mV)/mV)/second : Hz

m_inf = alpha_m/(alpha_m+beta_m) : 1

dh/dt = alpha_h*(1-h) - beta_h*h : 1

alpha_h = 32.4*exp(-0.14*(v+93.4*mV)/mV)/second : Hz
beta_h = 709/(1 + 4.2*exp(-0.06*(v+45.4*mV)/mV))/second : Hz

g_Na : siemens
'''

# eqs_i_Na = '''
# i_Na = (m**3) * h * g_Na * (v-E_Na) : amp
# 
# dm/dt = alpha_m * (1-m) - beta_m * m : 1
# 
# alpha_m = -0.1*(v+45*mV)/(exp((v+45*mV)/(-10*mV))-1)/(mV*second) : Hz
# beta_m = 4*exp((v+70*mV)/(18*mV))/second : Hz
# 
# dh/dt = alpha_h*(1-h) - beta_h*h : 1
# 
# alpha_h = 32.4*exp(-0.14*(v+93.4*mV)/mV)/second : Hz
# beta_h = 709/(1 + 4.2*exp(-0.06*(v+45.4*mV)/mV))/second : Hz
# 
# g_Na : siemens
# 
# '''

# Units on i_KK?
eqs_i_K = '''
i_K = i_KK + i_KNa : amp

i_KK = x * K_K * (K_o ** 0.59) * (K_i - K_o * exp(-1*v*F/(R*T))) : amp 

i_KNa = x * K_K * P_KNa * (K_o**0.59) * (Na_i - Na_o*exp(-1*v*F/(R*T))) : amp

dx/dt = (x_inf - x)/tau_x : 1
x_inf = 1/(1+exp((v+25.1*mV)/(-7.4*mV))) : 1
tau_x = 1*second/(17*exp(0.0398*v/mV) + 0.211*exp(-0.051*v/mV)) : second

P_KNa : 1
'''

eqs_i_f = '''
i_f = i_fNa + i_fK : amp

i_fNa = y*((K_o**1.83)/((K_o**1.83) + (K_mf**1.83)))*(g_fNa*(v-E_Na)) : amp
i_fK  = y*((K_o**1.83)/((K_o**1.83) + (K_mf**1.83)))*(g_fK*(v-E_K)) : amp

dy/dt = alpha_y * (1-y) - beta_y*y : 1
alpha_y = 0.36*(v+137.8*mV)/(exp(0.066*(v + 137.8*mV)/mV) - 1)/(mV*second) : Hz
beta_y  = 0.1*(v+76.3*mV)/(1-exp(-0.21*(v+76.3*mV)/mV))/(mV*second) : Hz

g_fNa : siemens
g_fK : siemens
K_K : amp*mM**-1.59
K_mf : mM
'''

eqs_i_p = '''
i_p = i_pmax*(Na_i/(Na_i + K_mNa))*(K_o/(K_o + K_mK))*(1-((v - 40*mV)/(211*mV))**2) : amp

i_pmax : amp
K_mK : mM
K_mNa : mM
'''

# Check units for i_NaCa
eqs_i_NaCa = '''
i_NaCa = k_NaCa*(x_2*k_21 - x_1*k_12)/(x_1 + x_2 + x_3 + x_4) : amp

x_1 = k_41*k_34*(k_23+k_21) + k_21*k_32*(k_43+k_41) : 1
x_2 = k_32*k_43*(k_14+k_12) + k_41*k_12*(k_34+k_32) : 1
x_3 = k_14*k_43*(k_23+k_21) + k_12*k_23*(k_43+k_41) : 1
x_4 = k_23*k_34*(k_14+k_12) + k_14*k_21*(k_34+k_32) : 1

k_43 = Na_i / (K_3ni + Na_i) : 1
k_12 = (Ca_i/K_ci)*exp(-1*Q_ci*v*F/(R*T))/d_i : 1
k_14 = ((Na_i**2)/(K_1ni*K_2ni)+(Na_i**3)/(K_1ni*K_2ni*K_3ni))*exp(Q_n*v*F/(2*R*T))/d_i : 1
k_41 = exp(-1*Q_n*v*F/(2*R*T)) : 1

d_i = 1 + Ca_i/K_ci + (Ca_i/K_ci)*exp(-1*Q_ci*v*F/(R*T)) + (Ca_i*Na_i)/(K_ci*K_cni) + Na_i/K_1ni + (Na_i**2)/(K_1ni*K_2ni) + (Na_i**3)/(K_1ni*K_2ni*K_3ni) : 1

k_34 = Na_o/(K_3no + Na_o) : 1
k_21 = (Ca_o/K_co)*exp(Q_co*v*F/(R*T))/d_o : 1
k_23 = ((Na_o**2)/(K_1no*K_2no) + (Na_o**3)/(K_1no*K_2no*K_3no))*exp(-1*Q_n*v*F/(2*R*T))/d_o : 1
k_32 = exp(Q_n*v*F/(2*R*T)) : 1

d_o = 1 + Ca_o/K_co + (Ca_o/K_co)*exp(Q_co*v*F/(R*T)) + Na_o/K_1no + (Na_o**2)/(K_1no*K_2no) + (Na_o**3)/(K_1no*K_2no*K_3no) : 1

K_1ni : mM
K_1no : mM
K_2ni : mM
K_2no : mM
K_3ni : mM
K_3no : mM
K_ci : mM
K_co : mM
K_cni : mM
k_NaCa : amp
Q_ci : 1
Q_co : 1
Q_n : 1
R : joule/mole/kelvin
T : kelvin
'''

eqs_i_bNa = '''
i_bNa = g_bNa*(v-E_Na) : amp

g_bNa : siemens
'''

eqs_i_bK = '''
i_bK = K_bK*(K_o**0.41)*(K_i - K_o*exp(-1*v*F/(R*T))) : amp

K_bK : amp*mM**-1.41
'''

eqs_i_up = '''
i_up = i_upmax*((Ca_i**2)/((Ca_i**2)+(K_mCaup**2))) : amp
i_tr = alpha_tr*Ca_up : amp
i_rel = alpha_rel*Ca_rel*((Ca_i**2)/((Ca_i**2) + (K_mCarel**2))) : amp

alpha_rel = 2*F*V_rel/tau_rel : amp/mM
alpha_tr  = 2*F*V_rel/tau_tr : amp/mM

i_upmax : amp
tau_rel : second
tau_tr  : second
K_mCarel : mM
K_mCaup : mM
'''

eqs_ion = '''
dNa_i/dt   = -1*(i_bNa + i_fNa + i_Na + 3*i_p + 3*i_NaCa + i_KNa)/(F*V_i) : mM
dNa_o/dt   = ((i_bNa + i_fNa + i_Na + 3*i_p + 3*i_NaCa + i_KNa)/(F*V_e)) + ((Na_b - Na_o)/tau_b) : mM

dK_i/dt    = -1*(i_KK + i_fK - 2*i_p)/(F*V_i) : mM
dK_o/dt    = ((i_KK + i_fK - 2*i_p)/(F*V_e)) + (K_b - K_o)/tau_b : mM

dCa_i/dt   = -1*(i_CaL + i_CaT - 2*i_NaCa + i_up - i_rel)/(2*F*V_i) : mM
dCa_o/dt   = (i_CaL + i_CaT - 2*i_NaCa)/(2*F*V_e) + (Ca_b - Ca_o)/tau_b : mM

dCa_up/dt  = (i_up-i_tr)/(2*F*V_up) : mM
dCa_rel/dt = (i_tr-i_rel)/(2*F*V_rel) : mM

tau_b : second

Ca_b : mM
K_b : mM
Na_b : mM
F : coulomb/mole
V_e : meter**3
V_i : meter**3
V_rel : meter**3
V_up : meter**3
'''

eqs_nernst = '''
E_Na = (R*T/(1*F))*log(Na_o/Na_i) : volt
E_K  = (R*T/(1*F))*log(K_o/K_i) : volt
E_Ca  = (R*T/(2*F))*log(Ca_o/Ca_i) : volt
'''

eqs = '''
dv/dt = -(i_tot)/Cm : volt
i_tot = i_CaL + i_CaT + i_Na + i_K + i_f + i_p + i_NaCa + i_bNa + i_bK + I_in : amp

I_in : amp
Cm : farad
'''

eqs += (eqs_i_CaL + eqs_i_CaT + eqs_i_Na + eqs_i_K + eqs_i_f + eqs_i_p + eqs_i_NaCa + eqs_i_bNa + eqs_i_bK + eqs_i_up + eqs_ion + eqs_nernst)

# Set Group
G = NeuronGroup(num_neurons, eqs,
                    threshold='v > 0*mV',
                    refractory='v > -20*mV',
                    method='euler')
                    
# Set initial conditions and constant values
G.alpha_f_2L = 3/second
G.beta_f_2L  = 40000/second/mM

G.tau_b = 0.1*second
G.tau_rel = 0.005*second
G.tau_tr = 0.4*second

G.Ca_b = 2*mM
G.K_b = 5.4*mM
G.Na_b = 140*mM

G.Cm = 52*pF # ALTERED FROM 32pF
G.F = faraday_constant

G.g_bNa = 0.24*nS
G.g_CaL = 400*nS
G.g_CaT = 85*nS
G.g_fK = 13.5*nS
G.g_fNa = 8.1*nS
G.g_Na = 250*nS

G.i_pmax = 226*pA
G.i_upmax = 21.2*pA

G.K_1ni = 395.3*mM
G.K_1no = 1628*mM
G.K_2ni = 2.289*mM
G.K_2no = 561.4*mM
G.K_3ni = 26.44*mM
G.K_3no = 4.663*mM
G.K_bK = 0.07*pA*mM**-1.41
G.K_ci = 0.0207*mM
G.K_cni = 26.44*mM
G.K_co = 3.663*mM
G.K_K = 0.26*pA*mM**-1.59
G.K_mCarel = 0.001*mM
G.K_mCaup = 0.0005*mM
G.K_mf = 10.3*mM
G.K_mK = 1*mM
G.K_mNa = 40*mM

G.k_NaCa = 4000*pA

G.P_KNa = 0.035

G.Q_ci = 0.1369
G.Q_co = 0
G.Q_n = 0.4315

G.R = gas_constant
G.T = 310*kelvin

G.V_e = 0.5*pliter
G.V_i = 2.5*pliter
G.V_rel = 0.015*pliter
G.V_up = 0.035*pliter
G.v = -64.9*mV

G.Ca_i = 0.000034*mM
G.Ca_up = 0.5832*mM
G.Ca_rel = 0.1101*mM
G.Ca_o = 2.0004*mM
G.K_i = 140.0073*mM
G.K_o = 5.4243*mM
G.Na_i = 7.4994*mM
G.Na_o = 139.9929*mM

G.y = 0.0287
G.x = 0.5682
G.m = 0.0139
G.h = 0.0087
G.d_L = 0.0001
G.f_L = 0.1505
G.f_2L = 0.2190
G.d_T = 0.0010
G.f_T = 0.1328

G.I_in = -00.0*pA

spikes         = SpikeMonitor(G)
currents       = StateMonitor(G,('v','i_CaL','i_CaT','i_Na','i_K','i_f','i_p','i_NaCa','i_bNa','i_bK','i_tot'),record=True)
# concentrations = StateMonitor(G,('Na_i','Na_o','K_i','K_o','Ca_i','Ca_o','Ca_up','Ca_rel'),record=True)
# #monitorCaL     = StateMonitor(G,('d_L','f_L','f_2L'),record=True)
# nernst         = StateMonitor(G,('E_Na','E_K','E_Ca'),record=True)
# i_Na_var       = StateMonitor(G,('alpha_m','beta_m','m','h','m_inf'),record=True)
# i_CaL_var      = StateMonitor(G,('d_L','f_L','f_2L'),record=True)

run(duration, report='text')


figure(1)
plot(currents.t/ms, currents.v[0]/mV)
xlabel('t (ms)')
ylabel('V (mV)')

figure(2)
plot(currents.t/ms, currents.i_tot[0]/pA)
xlabel('t (ms)')
ylabel('A (pA)')

# #figure(3)
# #plot(monitorCaL.t/ms, monitorCaL.d_L[0]/mM)
# #plot(monitorCaL.t/ms, monitorCaL.f_L[0]/mM)
# #plot(monitorCaL.t/ms, monitorCaL.f_2L[0]/mM)
# #xlabel('t (ms)')
# #ylabel('A (pA)')
# 
# figure(4)
# plot(concentrations.t/ms, concentrations.K_o[0]/mM)
# xlabel('t (ms)')
# ylabel('Conc. (mM)')
# ylim(0,200)
# 
# figure(5)
# plot(nernst.t/ms, nernst.E_Na[0]/mV)
# plot(nernst.t/ms, nernst.E_Ca[0]/mV)
# plot(nernst.t/ms, nernst.E_K[0]/mV)
# xlabel('t (ms)')
# ylabel('Pot. (mV)')
# 
# # figure(3)
# # plot(monitor2.t/ms, monitor2.K_o[0]/mM)
# # xlabel('t (ms)')
# # ylabel('Concentration (mM)')
# # ylim(0,6)
# 
# figure(6)
# #plot(i_Na_var.t/ms, i_Na_var.alpha_m[0])
# plot(i_Na_var.t/ms, i_Na_var.beta_m[0])
# #plot(i_Na_var.t/ms, i_Na_var.h[0])
# #plot(nernst.t/ms, nernst.K[0]/mV)
# xlabel('t (ms)')
# 
# figure(7)
# plot(i_Na_var.t/ms, i_Na_var.h[0])
# xlabel('t (ms)')
# 
# figure(8)
# plot(i_CaL_var.t/ms, i_CaL_var.d_L[0])
# #plot(i_CaL_var.t/ms, i_CaL_var.f_L[0])
# #plot(i_CaL_var.t/ms, i_CaL_var.f_2L[0])
# xlabel('t (ms)')
# ylabel('gating param (0-1)')

show()