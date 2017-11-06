# Will Huffman
# BME 503 Final Project
# Model based on Tao 2011 paper
# "A Model of Cellular Cardiac-Neural Coupling That Captures the Sympathetic Control
# of Sinoatrial Node Excitability in Normotensive and Hypertensive Rats"

# Altered to Parasympathetic Control Mechanisms

from brian2 import *

num_neurons = 1
duration = 2*second

# Parameters
area = 20000*umetre**2
Cm = 1*ufarad*cm**-2 * area
rest_pot = -60.0*mV
El = 10.613*mV+rest_pot
ENa = 115.0*mV+rest_pot
EK = -12.0*mV+rest_pot

gl = 0.3*msiemens*cm**-2 * area
gNa = 120*msiemens*cm**-2 * area
gK = 36*msiemens*cm**-2 * area

eqs_i_CaL = '''
i_CaL = d_L * f_L * f_2L * g_CaL * (v - E_Ca + 75) : amp

d_Linf = 1/(1+exp((v+6.6*mV)/(-6.6*mV))) : 1
tau_d_L = 0.002 * ms : second
dd_L/dt   = (d_Linf - d_L)/tau_d_L : 1

f_Linf = 1/(1+exp((v+25*mV)/(6*mV))) : 1
tau_f = 0.031*ms + 1*ms/(1+exp((v+37.6*mV)/(8.1*mV))) : second
df_L/dt = (f_Linf - f_L)/tau_f_L : 1

df_2L/dt = alpha_f_2L * (1 - f_2L) - beta_f_2L*Ca_i*f_2L : 1

alpha_f_2L : Hz
beta_f_2L : Hz
'''

eqs_i_CaT = '''
i_CaT = d_T * f_T * g_CaT * (v-E_CaT) : amp

d_T_inf = 1/(1+exp((v + 23*mV)/(-6.1*mV))) : 1
tau_d_T = 0.0006*ms + 0.0054*ms/(1+exp(0.03*(v + 100*mV)/mV)) : second
dd_T/dt = (d_Tinf - d)/tau_d_T : 1

f_Tinf = 1/(1+exp((v+75*mV)/(6.6*mV))) : 1
tau_f_T = 0.001*ms + 0.04*ms/(1+exp(0.08*(v+65*mV)/mV)) : second
df_T/dt = (f_Tinf - f_T)/tau_f_T : 1
'''

eqs_i_Na = '''
i_Na = (m**3) * h * g_Na * (v-E_Na) : amp

dm/dt = alpha_m * (1-m) - beta_m * m : 1
alpha_m = 200*(v + 34.3*mV)/(1-exp(-0.09*(v+34.3*mV)/mV))/ms : Hz 
beta_m = 8000*exp(-0.15*(v+56.2*mV)/mV)/ms : Hz

dh/dt = alpha_h*(1-h) - beta_h*h : 1
alpha_h = 32.4*exp(-0.14*(v+93.4*mV)/mV)/ms : Hz
beta_h = 709/(1 + 4.2*exp(-0.06*(v+45.4*mV)/mV))/ms : Hz
'''

# Units on i_KK?
eqs_i_K = '''
i_K = i_KK + i_KNa : amp

i_KK = x * K_K * (K_o ** 0.59) * (K_i - K_o * exp(-1*v*F/(R*T))) : amp 

i_KNa = x*K_K * P_KNa * (K_o**0.59)*(Na_i - Na_o*exp(-1*v*F/(R*T))) : amp

dx/dt = (x_inf - x)/tau_x : 1
x_inf = 1/(1+exp((v+25.1*mV)/(-7.4*mV))) : 1
tau_x = 1/(17*exp(0.0398*v/mV) + 0.211*exp(-0.051*v/mV))/ms : Hz
'''

eqs_i_f = '''
i_f = i_fNa + i_fK : amp

i_fNa = y*((K_o**1.83)/((K_o**1.83) + (K_mf**1.83)))*(g_fNa*(v-E_Na)) : amp
i_fK  = y*((K_o**1.83)/((K_o**1.83) + (K_mf**1.83)))*(g_fK*(v-E_K)) : amp

dy/dt = alpha_y * (1-y) - beta_y*y : 1
alpha_y = 0.36*(v+137.8*mV)/(exp(0.066*(v + 137.8*mV)/mV) - 1)/(mV*ms) : Hz
beta_y  = 0.1*(v+76.3*mV)/(1-exp(-0.21*(v+76.3*mV)/mV))/(mV*ms) : Hz
'''

eqs_i_p = '''
i_p = i_pmax*(Na_i/(Na_i + K_mNa))*(K_o/(K_o + K_mK))*(1-((v - 40*mV)/(211*mV))**2) : amp
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

d_i = 1 + Ca_i/K_ci + Ca_i/K_ci*exp(-1*Q_ci*v*F/(R*T)) + Ca_i*Na_i/(K_ci*K_cni) + Na_i/K_1ni + (Na_i**2)/(K_1ni*K_2ni) + (Na_i**3)/(K_1ni*K_2ni*K_3ni) : 1

k_34 = Na_o/(K_3no + Na_o) : 1
k_21 = (Ca_o/K_co)*exp(Q_co*v*F/(R*T))/d_o : 1
k_23 = ((Na_o**2)/(K_1no*K_2no) + (Na_o**3)/(K_1no*K_2no*K_3no))*exp(-1*Q_n*v*F/(2*R*T))/d_o : 1
k_32 = exp(Q_n*v*F/(2*R*T)) : 1

d_o = 1 + Ca_o/K_co + (Ca_o/K_co)*exp(Q_co*v*F/(R*T)) + Na_o/K_1no + (Na_o**2)/(K_1no*K_2no) + (Na_o**3)/(K_1no*K_2no*K_3no) : 1
'''

eqs_i_bNa = '''
i_bNa = g_bNa*(v-E_Na) : amp
'''

# Check units on i_bK
eqs_i_bK = '''
i_bK = K_bK*(K_o**0.41)*(K_i - K_o*exp(-1*v*F/(R*T))) : amp
'''

# Check units on all i equations
eqs_i_up = '''
i_up = i_upmax*((Ca_i**2)/((Ca_i**2)+(K_mCaup**2))) : amp
i_tr = alpha_tr*Ca_up : amp
i_rel = alpha_rel*Ca_rel*((Ca_i**2)/((Ca_i**2) + (K_mCarel**2))) : amp

alpha_rel = 2*F*V_rel/tau_rel : 1
alpha_tr  = 2*F*V_rel/tau_tr : 1
'''

eqs_ion = '''
dNa_i/dt   = -1*(i_bNa+i_fNa+i_Na+3*i_p + 3*i_NaCa + i_KNa)/(F*V_i) : 1
dNa_o/dt   = (i_bNa+i_fNa+i_Na+3*i_p + 3*i_NaCa + i_KNa)/(F*V_e) + (Na_b - Na_o)/tau_b : 1

dK_i/dt    = -1*(i_KK+i_fK-2*i_p+i_KAch)/(F*V_i) : 1
dK_o/dt    = (i_KK+i_fK-2*i_p+i_KAch)/(F*V_e) + (K_b - K_o)/tau_b : 1

dCa_i/dt   = -1*(i_CaL+i_CaT-2*i_NaCa+i_up-i_rel)/(2*F*V_i) : 1
dCa_o/dt   = (i_CaL+i_CaT-2*i_NaCa+i_up-i_rel)/(2*F*V_e) + (Ca_b - Ca_o)/tau_b : 1

dCa_up/dt  = (i_up-i_tr)/(2*F*V_up) : 1
dCa_rel/dt = (i_tr-i_rel)/(2*F*V_rel) : 1
'''

eqs = '''
dv/dt = -(i_tot)/Cm : volt
i_tot = i_CaL + i_CaT + i_Na + i_K + i_f + i_p + i_NaCa + i_bNa + i_bK : amp
'''

eqs += (eqs_i_CaL + eqs_i_CaT + eqs_i_Na + eqs_i_K + eqs_i_f + eqs_i_p + eqs_i_NaCa + eqs_i_bNa + eqs_i_bK + eqs_i_up + eqs_ion)

# Set Group
G = NeuronGroup(num_neurons, eqs,
                    threshold='v > -20*mV',
                    refractory='v > -20*mV',
                    method='exponential_euler')
                    
# Set initial conditions and constant values
G.alpha_f_2L = 3/second
G.beta_f_2L  = 40000/second


monitor = SpikeMonitor(G)
monitor2= StateMonitor(G,'v',record=True)

run(duration, report='text')


figure(1)
plot(monitor2.t/ms, monitor2.v[0]/mV) #plot the voltage for neuron 0 (index starts at 0)
xlim(0,25)
xlabel('t (ms)')
ylabel('V (mV)')

show()