####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Okinawa Institute of Science and Technology, Japan.
#
# This script runs on STEPS 2.x http://steps.sourceforge.net
#
# H Anwar, I Hepburn, H Nedelescu, W Chen and E De Schutter
# Stochastic calcium mechanisms cause dendritic calcium spike variability
# J Neuroscience 2013
#
# constants_withampa.py : provides a set of parameters and other constants 
# for the synaptically-induced dendritic ca burst model in the above study. 
# It is intended that this file is not altered.
#
# Script authors: Haroon Anwar and Iain Hepburn
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import math

# # # # # # # # # # # # # # # # SIMULATION CONTROLS # # # # # # # # # # # # # 

EF_DT = 5.0e-6          # The EField dt
NTIMEPOINTS =  3001

TIMECONVERTER =  2.0e-5

NITER = 1

############################PARAMETERS################################

init_pot = -60e-3

TEMPERATURE = 34.0

Q10 = 3

# # Faraday constant: unit of FARADAY is C/mol
# # Source: http://physics.nist.gov/cgi-bin/cuu/Value?f 24/2/2012
FARADAY = 96485.3365
#
# # Molar Gas Constant: unit of R is J/mol K
# # Source: http://physics.nist.gov/cgi-bin/cuu/Value?r 24/2/2012
R = 8.3144621
#
# # Avogadro constant: unit of AVOGADRO is /mol
# # Source: http://physics.nist.gov/cgi-bin/cuu/Value?na 24/2/2012
AVOGADRO = 6.02214129e23
#
# # Elementary charge: unit of E_CHARGE is C
# # Source: http://physics.nist.gov/cgi-bin/cuu/Value?e 24/2/2012
E_CHARGE = 1.602176565e-19


#FOR MSLO, THERE IS A NEW VALUE FOR Qt wrt to 25 degC
Qt = math.pow(Q10, ((TEMPERATURE-23)/10))
Qt_mslo = math.pow(Q10, ((TEMPERATURE-25)/10))

########## BULK RESISTIVITY ##########

Ra = 235.7*1.0e-2*100 # MODIFICTIONS FOR SMALL MODEL

########## MEMBRANE CAPACITANCE ##########

Cm = 0.64e-2
memb_capac_proximal = 1.2*Cm
memb_capac_spiny = 5.3*Cm


########## CaP channels density & permiability per channel ##########

# CaP_P is permiability per channel (m3/s)
# CaP_ro is channel/surface area (/m2)
# P in Ca Dynamics model is 0.95e-4 cm/s --> 0.95e-6 m/s

CaP_P = 2.5e-20 
CaP_ro = 3.8e13

##########CaP channel parameters ####################

#Units (mV)
vhalfm = -29.458
cvm = 8.429

def minf_cap(V):
    #Units (mV)
    vhalfm = -29.458
    cvm = 8.429
    vshift = 0.0
    
    return (1.0/(1.0 + math.exp(-(V-vhalfm-vshift)/cvm)))

def tau_cap(V):
    vshift = 0.0
    if (V-vshift) >= -40:
        return (0.2702 + 1.1622 * math.exp(-(V+26.798-vshift)*(V+26.798-vshift)/164.19))
    else:
        return (0.6923 * math.exp((V-vshift)/1089.372))

def alpha_cap(V):
    return (minf_cap(1e3 * V)/tau_cap(1e3 * V)) * Qt * 1e3

def beta_cap(V):
    return ((1.0-minf_cap(V * 1e3))/tau_cap(V * 1e3)) * Qt * 1e3


## Intitial conditions

CaP_m0_p = 0.92402
CaP_m1_p = 0.073988
CaP_m2_p = 0.0019748
CaP_m3_p = 1.7569e-05


########## CaT channels density & permiability per channel ##########

# CaT_P is permiability per channel (m3/s)
# CaT_ro is channel/surface area (/m2)
# P in Ca Dynamics model is 6.2e-6 cm/s -->6.2e-8 m/s
# P in Ca Dynamics model with ampa is 3.24e-6 cm/s --> 3.24e-8 m/s

CaT_P = 1.65e-20
CaT_ro = 1.9636e12

#CaT_ro = 3.7576e12 (previously used value in model with no ampa)


def minf_cat(V):
    #Units (mV)
    vhalfm = -52.0
    cvm = -5.0
    vshift = 0.0
    
    return (1.0/(1.0 + math.exp((V-vhalfm-vshift)/cvm)))

def taum_cat(V):
    vshift = 0.0
    if V > -90.0:
        return (1.0 + 1.0 / (math.exp((V+40.0-vshift)/9.0) + math.exp(-(V+102.0-vshift)/18.0)))
    else:
        return 1.0

def hinf_cat(V):
    vhalfh = -72.0
    cvh = 7.0
    vshift = 0.0
    return (1.0/(1.0 + math.exp((V-vhalfh-vshift)/cvh)))

def tauh_cat(V):
    vshift = 0.0
    return (15.0 + 1.0 / (math.exp((V+32.0-vshift)/7.0)))

def alpham_cat(V):
    return (minf_cat(V)/taum_cat(V))

def betam_cat(V):
    return ((1-minf_cat(V))/taum_cat(V))

def alphah_cat(V):
    return (hinf_cat(V)/tauh_cat(V))

def betah_cat(V):
    return ((1-hinf_cat(V))/tauh_cat(V))

## Initial conditions

CaT_m0h0_p = 0.58661
CaT_m1h0_p = 0.23687
CaT_m2h0_p = 0.023912
CaT_m0h1_p = 0.10564
CaT_m1h1_p = 0.042658
CaT_m2h1_p = 0.0043063

########## BK channels density & conductance per channel ##########

# Total conductance = BK_G (conductance/channel) * BK_ro (channel/surface area)
# BK in Ca Dynamics model is 4.25e-2 S/cm2 --> 4.25e2 S/m2

BK_G = 2.1e-10
BK_ro = 2.0238e12 * 1.5
BK_rev = -77e-3

######### BK channel parameters ######################

#Units (1)
Qo = 0.73
Qc = -0.67

#Units (/s)
pf0 = 2.39
pf1 = 5.4918
pf2 = 24.6205
pf3 = 142.4546
pf4 = 211.0220

pb0 = 3936
pb1 = 687.3251
pb2 = 234.5875
pb3 = 103.2204
pb4 = 11.6581

#Units(/M)
k1 = 1.0e6

#Units(/s)
onoffrate = 1.0e3

L0 = 1806

#Units (M)
Kc = 8.63e-6
Ko = 0.6563e-6


c_01 = 4.*k1*onoffrate*Qt_mslo
c_12 = 3.*k1*onoffrate*Qt_mslo
c_23 = 2.*k1*onoffrate*Qt_mslo
c_34 = 1.*k1*onoffrate*Qt_mslo
o_01 = 4.*k1*onoffrate*Qt_mslo
o_12 = 3.*k1*onoffrate*Qt_mslo
o_23 = 2.*k1*onoffrate*Qt_mslo
o_34 = 1.*k1*onoffrate*Qt_mslo

c_10 = 1.*Kc*k1*onoffrate*Qt_mslo
c_21 = 2.*Kc*k1*onoffrate*Qt_mslo
c_32 = 3.*Kc*k1*onoffrate*Qt_mslo
c_43 = 4.*Kc*k1*onoffrate*Qt_mslo
o_10 = 1.*Ko*k1*onoffrate*Qt_mslo
o_21 = 2.*Ko*k1*onoffrate*Qt_mslo
o_32 = 3.*Ko*k1*onoffrate*Qt_mslo
o_43 = 4.*Ko*k1*onoffrate*Qt_mslo


f_0 = lambda mV: pf0*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
f_1 = lambda mV: pf1*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
f_2 = lambda mV: pf2*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
f_3 = lambda mV: pf3*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
f_4 = lambda mV: pf4*Qt_mslo*(math.exp((Qo* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))

b_0 = lambda mV: pb0*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
b_1 = lambda mV: pb1*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
b_2 = lambda mV: pb2*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
b_3 = lambda mV: pb3*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))
b_4 = lambda mV: pb4*Qt_mslo*(math.exp((Qc* FARADAY* mV) / (R* (TEMPERATURE + 273.15))))


# Initial conditions
BK_C0_p= 0.99997
BK_C1_p= 4.3619e-07
BK_C2_p= 4.1713e-09
BK_C3_p= 4.4449e-11
BK_C4_p= 6.3132e-14

BK_O0_p= 2.5202e-05
BK_O1_p= 1.1765e-06
BK_O2_p= 6.6148e-08
BK_O3_p= 2.4392e-09
BK_O4_p= 4.0981e-11

########## SK channel density & conductance per channel #############

# Total conductance = SK_G (conductance/channel) * SK_ro (channel/surface area)
# SK in Ca Dynamics model is 3.1e-4 S/cm2 --> 3.1 S/m2

SK_G = 1.0e-11
SK_ro = 31.0e10

SK_rev = -77e-3

######### SK channel parameters ###################

#Units (/s)
invc1 = 80
invc2 = 80
invc3 = 200

invo1 = 1000
invo2 = 100

diro1 = 160
diro2 = 1200

#Units ( /s M)

dirc2 = 200e6
dirc3 = 160e6
dirc4 = 80e6

invc1_t = invc1*Qt
invc2_t = invc2*Qt
invc3_t = invc3*Qt

invo1_t = invo1*Qt
invo2_t = invo2*Qt

diro1_t = diro1*Qt
diro2_t = diro2*Qt

dirc2_t = dirc2*Qt/3.0
dirc3_t = dirc3*Qt/3.0
dirc4_t = dirc4*Qt/3.0


# Intital conditions
SK_C1_p= 0.96256
SK_C2_p= 0.036096
SK_C3_p= 0.0010829
SK_C4_p= 6.4973e-06

SK_O1_p= 0.00017326
SK_O2_p= 7.7967e-05

######## AMPA rate constants ##############
#Total conductance = 20nS, 30nS and 40nS ---> 20e-9 S, 30e-9 S and 40e-9 S
#Single AMPA receptor conductance (Hausser and Roth 1997; Momiyama et al. 2003; Tanaka et al. 2005) - 7-8 pS

#Units (S)

AMPA_G = 7e-12
AMPA_TotalG = 500e-9*0.01 # MODIFIED FOR SMALL MODEL

#Units (1)

AMPA_receptors = AMPA_TotalG/AMPA_G

#Units (V)

AMPA_rev = 0e3

#Units (/s M)

rb = 13e6

#Units (/s)

ru1 = 0.0059e3
ru2 = 86e3
ro = 2.7e3
rc = 0.2e3
rd = 0.9e3
rr = 0.064e3


######### leak current channel density & conductance per channel ########
# Total conductance = 1e-6 S/cm2 --> 1e-2 S/m2
# Total membrane resistance  = 122 kohm.cm2 

Rm = 122e3/1.0e4
g_leak_proximal = 1.2/Rm
g_leak_spiny = 5.3/Rm


L_ro_proximal = 25.0e10

L_G =  g_leak_proximal/L_ro_proximal

L_ro_spiny = g_leak_spiny/L_G


L_rev = -61e-3


######### Pump parameters ###################

P_f_kcst = 3e9
P_b_kcst = 1.75e4
P_k_kcst = 7.255e4


############################CALCIUM BUFFERING MODEL################################

########## Ca concentrations #########

Ca_oconc = 2e-3
Ca_iconc = 45e-9

########## Mg concentrations #########

Mg_conc = 590e-6

########## Buffer concentrations #############

iCBsf_conc = 27.704e-6
iCBCaf_conc = 2.6372e-6
iCBsCa_conc= 1.5148e-6
iCBCaCa_conc= 0.14420e-6

CBsf_conc= 110.82e-6
CBCaf_conc= 10.549e-6
CBsCa_conc= 6.0595e-6
CBCaCa_conc= 0.57682e-6

PV_conc= 3.2066e-6
PVCa_conc= 16.252e-6
PVMg_conc= 60.541e-6

# Diffusion constant of Calcium
DCST = 0.223e-9
# Diffusion constant of Calbindin (CB)
DCB = 0.028e-9
# Diffusion constant of Parvalbumin (PV)
DPV = 0.043e-9

#iCBsf-fast
iCBsf1_f_kcst = 4.35e7
iCBsf1_b_kcst = 35.8

#iCBsCa
iCBsCa_f_kcst = 0.55e7
iCBsCa_b_kcst = 2.6

#iCBsf_slow
iCBsf2_f_kcst = 0.55e7
iCBsf2_b_kcst = 2.6

#iCBCaf
iCBCaf_f_kcst = 4.35e7
iCBCaf_b_kcst = 35.8

#CBsf-fast
CBsf1_f_kcst = 4.35e7
CBsf1_b_kcst = 35.8

#CBsCa
CBsCa_f_kcst = 0.55e7
CBsCa_b_kcst = 2.6

#CBsf_slow
CBsf2_f_kcst = 0.55e7
CBsf2_b_kcst = 2.6

#CBCaf
CBCaf_f_kcst = 4.35e7
CBCaf_b_kcst = 35.8

#PVca
PVca_f_kcst = 10.7e7
PVca_b_kcst = 0.95

#PVmg
PVmg_f_kcst = 0.8e6
PVmg_b_kcst = 25

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

