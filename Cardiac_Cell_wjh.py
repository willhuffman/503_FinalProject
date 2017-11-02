# Will Huffman
# BME 503 Final Project
# Model based on Tao 2011 paper
# "A Model of Cellular Cardiac-Neural Coupling That Captures the Sympathetic Control
# of Sinoatrial Node Excitability in Normotensive and Hypertensive Rats"

# Altered to Parasympathetic Control Mechanisms

from brian2 import *

num_neurons = 100
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

eqs_iCaL = '''
iCaL = f_PKA_L * gCaL * (v-ECaL) * d * f * fCa : amp

d_bar   = 1/(1+ exp(-(v+12.1*mV + v_SHIFT_L)/(6*mV))) : 1
tau_d   = 1/(alpha_d + alpha_d) : second
alpha_d = (-0.02839 * (v + 35*mV))/(exp(-(v+35*mV)/(2.5*mV)) - 1) - (0.0849*v)/(exp(-v/(4.808*mV)) - 1) : Hz
beta_d  = 0.0143 * (v - 5*mV)/(exp((v-5*mV)/mV) - 1) : Hz
dd/dt   = (d_bar - d)/tau_d : 1

f_bar = 1/(1+exp(v(+30*mV+v_SHIFT_L)/(5*mV))) : 1
tau_f = 55.3*ms + 257.1*exp(-((v+32.5*mV)/(13.9*mV))**2) : second
df/dt = (f_bar - f)/tau_f : 1

alpha_fCa = Km_fCa * beta_fCa : Hz
fCa_bar = Km_fCa / (Km_fCa + Ca_sub) : 1
tau_fCa = fCa_bar / alpha_fCa : second
dfCa =  (fCa_bar - fCa) / tau_fCa : 1
'''

eqs_iCaT = '''
iCaT = g_CaT * (v-ECaT) * d * f : amp

d_bar = 1/(1+exp(-(v+26.3*mV)/(6*mv))) : 1
tau_d = 1/(1.068*exp((v+26.3*mV)/(30*mV)) + 1.068*exp(-(v+26.3*mV)/(30*mV))) : second
dd/dt = (d_bar - d)/tau_d : 1

f_bar = 1/(1+exp((v+61.7*mV)/(5.6*mV))) ; 1
tau_f = 1/(0.0153*exp(-(v+61.7*mV)/(5.6*mV)) + 0.015*exp((v+61.7*mV)/(15.38*mV))) : second
df/dt = (f_bar - f)/tau_f : 1
'''

eqs_iKr = '''
iKr = g_Kr * (v - EK) * (0.6*paF + 0.4*paS)*piy : amp

pa_bar = 1/(1+exp(-(v+24*mV)/(5*mV))) : 1

tau_paS = 1/(0.0042*exp((v-9*mV)/(17*mV)) + 0.00015*exp(-(v-9*mV)/(22.6*mV))) : second
dpaS/dt = (paS_bar - paS)/tau_paS : 1

tau_paF = 1/(0.0372*exp((v-9*mV)/(15.9*mV)) + 0.00096*exp(-(v-9*mV)/(22.5*mV))) : second
dpaF/dt = (paF_bar - paF)/tau_paF : 1

piy_bar = 1/(1+exp((v+9.6*mV)/(10.1*mV))) : 1
tau_piy = 2*ms : second
dpiy/dt = (piy_bar - piy)/tau_piy : 1
'''
eqs_iKs = '''
iKs = f_Ks * g_Ks * (v-EKs) * n**2 : amp

n_bar = alpha_n/(alpha_n + beta_n) : 1
tau_n = 1/(alpha_n + beta_n) : second
alpha_n = 0.014/(1+exp(-(v - 40*mV + SHIFT_Ks)/(45*mV))) : Hz
beta_n = 0.001*exp(-(v+SHIFT_Ks)/(22*mV)) : Hz
dn/dt = (n_bar - n)/tau_n : 1
'''

eqs_ito = '''
ito = g_to*r*q*(v - EK) : amp

r_bar = 1/(1+exp(-(v-19.3*mV)/(15*mV))) : 1
tau_r = 14.405/(1.037*exp(0.09*(v+39.61*mV)/mV) + 0.369*exp(-0.12*(v+23.84*mV)/mV)) + 2.7535 : second
dr/dt = (r_bar - r)/tau_r : 1

q_bar = 1/(1+exp((v+49*mV)/(13*mV))) : 1
tau_q = 39.102/(0.57*exp(-0.08*(v+44*mV)/mV)+0.065*exp(0.1*(v+45.93*mV)/mV)) + 6.06 : second
dq/dt = (q_bar - q)/tau_q : 1

tau_Sslow = 3.7*exp(-((v+70.0*mV/30.0)/mV)**2)*ms + 0.035*ms : second
'''

eqs_isus = '''
isus = g_sus * r * (v- EK) : amp
r_bar = 1/(1+exp(-(v-19.3*mV)/(15*mV))) : 1
tau_r = 14.405*ms/(1.037*exp(0.09*(v+30.61*mV)/mV) + 0.369*exp(-0.12*(v + 23.84*mV)/mV)) + 2.7535*ms : second
dr/dt = (r_bar - r)/tau_r : 1
'''

eqs_if = '''
if = g_f * ( ifNa + ifK) : amp

ifNa = g_fNa * (v - ENa) * y**2 : amp
IfK = g_K * (v - EK) * y**2 : amp

y_bar = 1/(1+exp((v+100*mv-v_shift)/(13.5*mV))) : 1
tau_y = 0.7166529*ms/(exp(-(v+425.5*mV)/(45.302*mV))+exp((v-73.08*mV)/(19.231*mV))) : 1
dy/dt = (y_bar - y)/tau_y : 1
'''

eqs_ist = '''
ist = f_PKA_st*g_st*(v-Est)*qa*qi : amp

qa_bar = 1/(1+exp(-(v+57*mv)/(5*mv)) : 1
tau_qa = 1/(alpha_qa + beta_qa) : second
alpha_qa = 1/(0.15*exp(-1*v/(11*mV))+0.2*exp(-1*v/(700*mV))) : Hz
beta_qa = 1/(16*exp(v/(8*mV))+15*exp(v/(mV))) : Hz
dqa/dt = (qa_bar - qa)/tau_qa : 1

qi_bar = alpha_qi/(alpha_qi + beta_qi) : 1
tau_qi = alpha_qi/(alpha_qi + beta_qi) : second
alpha_qi = 1/(3100*exp(v/(13*mV)) + 700*exp(-1*v/(70*mV))) : Hz
beta_qi = 1/(95*exp(-1*v/(10*mV)) + 50*exp(-1*v/(700*mV))) + 0.000229/(1+exp(-1*v/(5*mV))) : Hz
dqi/dt = (qi_bar - qi)/tau_qi : 1
'''

eqs_ibNa = '''
ibNa = g_bNa*(v-ENa)
'''

eqs_iKAch = '''
iKAch = g_KAch * (Ki - Ko*exp(-1*v*F/(R*T))
'''

eqs_iNaK = '''
iNaK = iNaK_max * (1+ (Km_Kp/Ko)**1.2)**-1 * (1+(Km_Nap/Nai)**1.3)**-1*(1+exp(-1*(v_ENa+120)/30))**-1 : amp
'''

eqs_NaCa = '''
iNaCa = k_NaCa * ((Nai**3)*Cao*exp(0.03743*v*lambda_NaCa/mV)-Nao**3*Casub*exp(0.03743*v*(lambda_NaCa-1)))/(1.0+d_NaCa*(Nai**3*Cao+Nao**3*Casub)) : amp
'''

eqs_Cai = '''
J_Cadiff = (Casub - Cai)/tau_diffCa : 1
J_rel = P_rel * (Carel - Casub)*(Casub**2)/(Casub**2+K_rel**2) : 1
J_up = P_up * Cai/(Cai + K_up) : 1
J_tr = (Caup - Carel)/tau_tr
'''

eqs_ion = '''
dCai/dt = (J_Cadiff*v_sub-J_up*v_up-i_CaP)/v_i -(CM_tot*f_CMi_rate+TC_tot*f_TC_rate+TMC_tot*f_TMC_rate) : 1
dCasub/dt = -1*((iCaL + iCaT - 2.0*iNaCa)*Cm/(2*F)+J_rel*v_rel)/v_sub - (J_Cadiff + CM_tot*f_CMs_rate) : 1
'''

#eqs_ina = '''
#ina = gNa * (m**3) * h * (ENa-v) :  amp
#
#dm/dt = alpham * (1-m) - betam*m : 1
#alpham = (0.1/mV) * (-(v-rest_pot)+25*mV) / (exp((-(v-rest_pot)+25*mV) / (10*mV)) - 1)/ms : Hz
#betam = 4*exp(-(v-rest_pot)/(18*mV))/ms : Hz
#
#dh/dt = alphah * (1-h) - betah*h : 1
#alphah = 0.07 * exp((-(v-rest_pot)) / (20*mV))/ms : Hz
#betah = 1 / (exp((-(v-rest_pot)+30*mV) / (10*mV)) + 1)/ms : Hz 
#'''

#eqs_ik = '''
#ik = gK * (n**4) * (EK-v):amp
#dn/dt = alphan * (1-n) - betan * n : 1
#alphan = (0.01/mV) * (-(v-rest_pot)+10*mV) / (exp((-(v-rest_pot)+10*mV) / (10*mV)) - 1)/ms : Hz
#betan = 0.125*exp((-(v-rest_pot))/(80*mV))/ms : Hz
#'''

#eqs_il = '''
#il = gl * (El-v) :amp
#'''

#eqs = '''
#dv/dt = (ina+ik+il +I)/Cm:  volt
#I : amp
#'''

eqs = '''
dv/dt = -(iCaL + iCaT + iKr + iKs + ist + ito + isus + if + iKAch + iBNa + iNaK + iNaCa)/Cm : volt

'''

eqs += (eqs_ina+eqs_ik+eqs_il) 

# Threshold and refractoriness are only used for spike counting
group = NeuronGroup(num_neurons, eqs,
                    threshold='v > -20*mV',
                    refractory='v > -20*mV',
                    method='exponential_euler')
group.v = rest_pot
group.m=0.0529
group.n=0.3177
group.h=0.596


monitor = SpikeMonitor(group)
monitor2= StateMonitor(group,'v',record=True)

#group.I = 0*nA
#run(5.0*ms,report='text')
#group.I = 1.5*nA
#run(1.0*ms)
#group.I = 0*nA
#run(20.0*ms)

group.I = '(7.0*nA * i) / num_neurons'
run(2000.0*ms, report='text')

figure(1)
plot(group.I/nA, monitor.count/duration)
xlabel('I (nA)')
ylabel('Firing rate (sp/s)')
xlim(0,7)
#ylim(0,200)

figure(2)
plot(monitor2.t/ms, monitor2.v[0]/mV) #plot the voltage for neuron 0 (index starts at 0)
plot(monitor2.t/ms, monitor2.v[50]/mV) #plot the voltage for neuron 0 (index starts at 0)
plot(monitor2.t/ms, monitor2.v[99]/mV) #plot the voltage for neuron 0 (index starts at 0)
xlim(0,25)
#ylim(-80,60) #set axes limits
xlabel('t (ms)')
ylabel('V (mV)')

show()