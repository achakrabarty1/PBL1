# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 01:30:54 2021

@author: stper
"""

import Humidifier as hum
import splitter as sp
import Metabolism as met
import capillaries as cap
from functions import volumetric_to_molar
#from functions import getToBlood
import numpy as np
import matplotlib.pyplot as plt



# data constants
SAT_PCT_O2 = 0.1967
SAT_PCT_CO2 = 0.003
SAT_PCT_H2O = 0.062
ALV_PCT_O2 = 0.136
ALV_PCT_CO2 = 0.053
ALV_PCT_H2O = 0.062
DIST = .03 #distance from alveola to capillary

def getGradient(pBlood, pAlv):
    return (pBlood-pAlv)/DIST


NORMAL_ALV_AREA = 8.5 * (10^5)  #"true" alveolar area in cm^2 from Gehr et. al, 1978
ALV_AREA = 3.0 * (10^5)  #cm^2          MODIFY THIS ONE FOR DATA!
PO2_BLOOD = 104 #mm hg
PO2_ALV = 40 #mm hg
PO2_GRADIENT = getGradient(PO2_BLOOD, PO2_ALV)
DIFF_COEFF_O2 = 2.4*(10^-5) #cm^2/s

PCO2_BLOOD = 40
PCO2_ALV = 45
PCO2_GRADIENT = getGradient(PCO2_BLOOD, PCO2_ALV)
DIFF_COEFF_CO2 = 1.55*(10^-5) #cm^2/s

mw_gluc = 180.16 #g/mol
c_gluc_in = 0.9 #g/L
V_b = 5.0 #L/min
m_H2O_in = 2100 #mL/day
mw_H2O = 18.015 #g/mol

def calcToBlood(coeff,area, gradient):
    return (coeff)*(area)*(gradient)

def getNormalArea():
    return NORMAL_ALV_AREA
def getArea():
    return ALV_AREA

#ALVEOLI
v_toBlood_O2 = calcToBlood(DIFF_COEFF_O2, ALV_AREA, PO2_GRADIENT)
n_toBlood_O2 = volumetric_to_molar(v_toBlood_O2)
v_toBlood_CO2 = calcToBlood(DIFF_COEFF_CO2, ALV_AREA, PCO2_GRADIENT)
n_toBlood_CO2 = volumetric_to_molar(v_toBlood_CO2)
#n_sat_alveoli = sp.getSplitterAlv()
#n_fromBlood_CO2 = met.get_CO2_met()
#n_alv = n_sat_alveoli

#CAPILLARIES
ki = NORMAL_ALV_AREA/ALV_AREA   #add another constant here? Not sure how this would scale
                                        #as i'm not sure how much SA is lost with emphysema

MW_Hb = 64500    #g HB/mol
CO = 5   #5 L blood/min
Dens_Hb = 150 #g/L blood
HbOrel = 0.59  #ideal L O2/(100 g Hb * min)

v_O2_alv_cap= v_toBlood_O2
n_O2_alv_cap = ki*n_toBlood_O2
n_O2_cap_met = n_O2_alv_cap  #mol/min
n_O2_met_cap = 0

n_Hb_ideal = (1/(HbOrel*100)/MW_Hb)*v_O2_alv_cap
normal_n_Hb_EPO = Dens_Hb*CO/(100*MW_Hb)
n_Hb_EPO = (ki)*normal_n_Hb_EPO
n_Hb_cap_met = Dens_Hb*CO/MW_Hb+n_Hb_EPO
n_HbO_cap_met = n_O2_alv_cap/4
n_Hbf_cap_met = n_Hb_cap_met - n_HbO_cap_met

n_CO2 = .0225*5   #mol/min
n_CO2_alv_cap = n_toBlood_CO2
n_CO2_met_cap = n_O2_cap_met
n_CO2_cap_alv = n_CO2_alv_cap + n_CO2_met_cap

#METABOLISM
n_Hb_met_cap = n_Hb_cap_met
n_Hb_waste = n_Hb_cap_met + n_Hb_EPO - n_Hb_met_cap
n_gluc_ingest = V_b * c_gluc_in / mw_gluc #mol/min
n_H2O_ingest = m_H2O_in * 24 * 60 / mw_H2O

#determine limiting reagent
LR = min(n_O2_cap_met/6, n_gluc_ingest)
print(LR)

# calculations of consumed products
# because we found O2 to be the limiting reagent,
n_O2_con = n_O2_cap_met
n_gluc_con = n_O2_con / 6
n_HbO_con = n_O2_con / 4

# calculations of generated products
# because we found O2 to be the limiting reagent, and RQ = 1:
n_CO2_gen = n_O2_cap_met
n_H2O_gen = n_O2_cap_met
n_Hbf_gen = n_HbO_con

# calculations of output
n_CO2_met_cap = n_CO2_gen
n_H2O_met_cap = n_H2O_ingest + n_H2O_gen
n_O2_met_cap = n_O2_cap_met - n_O2_con
n_gluc_met_cap = n_gluc_ingest - n_gluc_con
n_Hbf_met_cap = n_Hbf_cap_met + n_Hbf_gen
n_HbO_met_cap = n_HbO_cap_met - n_HbO_con



def getNBloodO2():
    return n_toBlood_O2
def getVBloodO2():
    return v_toBlood_O2
def getNBloodCO2():
    return n_toBlood_CO2
def getVBloodCO2():
    return v_toBlood_CO2
 
            

    
a0= NORMAL_ALV_AREA     #initial area
aEnd=  3.0 * (10^5)    #area of a severely emphysematic lung in cm^2 from Thurlbeck, 1967
da=-5     #small decrement in area
a=np.arange(a0,aEnd, da) #creates an array with area values from healthy to diseased area
                        #separated by da   
print(a)                        

n=int((aEnd-a0)/da)  #the number of steps it takes to get from healthy to diseased
print(n)
l = n+1
k = np.zeros(l)
vtbO=np.zeros(l)
ntbO=np.zeros(l)
vtbC=np.zeros(l)
ntbC=np.zeros(l)
nfbC=np.zeros(l)
k[0] = ki
vtbO[0]=v_toBlood_O2
ntbO[0]=n_toBlood_O2
vtbC[0]=v_toBlood_CO2
ntbC[0]=n_toBlood_CO2
#nfbC[0]=met.get_CO2_met()

vOac= np.zeros(l)
nOac = np.zeros(l)
nOcm= np.zeros(l)  #mol/min
nOmc = np.zeros(l)
nOc = np.zeros(l)
nCac = np.zeros(l)
nCmc = np.zeros(l)
nCac = np.zeros(l)
nCac = np.zeros(l)
nCca = np.zeros(l)
'''
vOac[0] = v_O2_alv_cap
nOac[0] = n_O2_alv_cap
nOcm[0] = n_O2_cap_met
nOmc[0] = n_O2_met_cap
nCac[0] = n_CO2_alv_cap
nCmc[0] = n_CO2_met_cap
nOc[0] = n_O2_con
nCca[0] = n_CO2_alv_cap + n_CO2_met_cap
'''

nHE = np.zeros(l)
nHcm = np.zeros(l)
nHOcm = np.zeros(l)
nHfcm = np.zeros(l)
'''
nHE[0] = n_Hb_EPO
nHcm[0] = n_Hb_cap_met
nHOcm[0] = n_HbO_cap_met
nHfcm[0] = n_Hbf_cap_met    
'''   
    
for i in range(0,len(a)-1):
    k[i] = NORMAL_ALV_AREA/a[i] 
    vtbO[i] = calcToBlood(DIFF_COEFF_O2, -a[i], getGradient(PO2_BLOOD, PO2_ALV))
    ntbO[i]=volumetric_to_molar(vtbO[i])
    vtbC[i]=calcToBlood(DIFF_COEFF_CO2, a[i], getGradient(PCO2_BLOOD, PCO2_ALV))
    ntbC[i]=volumetric_to_molar(vtbC[i])
    vOac[i] = vtbO[i]
    nOac[i] = ntbO[i]
    nOcm[i] = nOac[i]
    nOc[i] = nOcm[i]
    
    nOmc[i] = nOac[i] - nOc[i]
    nCac[i] = ntbC[i]
    nCmc[i] = nOcm[i]
    nCca[i] = ntbC[i] + nCmc[i]
    
    nHE[i+1] = k[i]*normal_n_Hb_EPO
    nHcm[i+1] = Dens_Hb*CO/MW_Hb+nHE[i]
    nHOcm[i+1] = nOac[i+1]/4
    nHfcm[i+1] = nHcm[i+1]-nHOcm[i+1]


n_CO2 = .0225*5   #mol/min

n_CO2_cap_alv = n_CO2_alv_cap + n_CO2_met_cap

#METABOLISM
n_Hb_met_cap = n_Hb_cap_met
n_Hb_waste = n_Hb_cap_met + n_Hb_EPO - n_Hb_met_cap
n_gluc_ingest = V_b * c_gluc_in / mw_gluc #mol/min
n_H2O_ingest = m_H2O_in * 24 * 60 / mw_H2O

#determine limiting reagent
LR = min(n_O2_cap_met/6, n_gluc_ingest)
print(LR)

# calculations of consumed products
# because we found O2 to be the limiting reagent,
n_O2_con = n_O2_cap_met
n_gluc_con = n_O2_con / 6
n_HbO_con = n_O2_con / 4

# calculations of generated products
# because we found O2 to be the limiting reagent, and RQ = 1:
n_CO2_gen = n_O2_cap_met
n_H2O_gen = n_O2_cap_met
n_Hbf_gen = n_HbO_con

# calculations of output
n_CO2_met_cap = n_CO2_gen
n_H2O_met_cap = n_H2O_ingest + n_H2O_gen
n_O2_met_cap = n_O2_cap_met - n_O2_con
n_gluc_met_cap = n_gluc_ingest - n_gluc_con
n_Hbf_met_cap = n_Hbf_cap_met + n_Hbf_gen
n_HbO_met_cap = n_HbO_cap_met - n_HbO_con
        
figure=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis=figure.add_subplot(1,1,1)
axis.grid(True)    
axis.plot(a, vtbO, label='volumetric O2 flow from Alv to Cap')
axis.legend()        
axis.set(xlabel='Area (cm^2)', ylabel='Oxygen flow (L/min)') 
figure.savefig('vtbo')  
        
figure2=plt.figure(num=2, clear=True) #plotting for arteriosclerosis v treated
axis2=figure2.add_subplot(1,1,1) 
axis2.grid(True)     
axis2.plot(a, ntbO, label='Molar O2 flow from Alv to Cap')       
axis2.legend()
axis2.set(xlabel='Area (cm^2)', ylabel='O2 flow (L/min)') 
figure2.savefig('ntbo')  

figure3=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis3=figure3.add_subplot(1,1,1)
axis3.grid(True)    
axis3.plot(a, vtbC, label='volumetric CO2 flow from Alv to Cap')
axis3.legend()        
axis3.set(xlabel='Area (cm^2)', ylabel='CO2 flow (L/min)') 
figure3.savefig('vtbc')

figure4=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis4=figure4.add_subplot(1,1,1)
axis4.grid(True)    
axis4.plot(a, ntbC, label='Molar CO2 flow from Alv to Cap')
axis4.legend()        
axis4.set(xlabel='Area (cm^2)', ylabel='CO2 flow (L/min)') 
figure4.savefig('ntbc')

figure5=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k, label='Hb coefficient')
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel='Hb coeff') 
figure5.savefig('k')

figure6=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis6=figure.add_subplot(1,1,1)
axis6.grid(True)    
axis6.plot(a, nOmc, label='Molar flow of Oxygen from Metabolic Reactor to Capillaries')
axis6.legend()        
axis6.set(xlabel='Area (cm^2)', ylabel='O2 (mol/min)') 
figure6.savefig('nOmc')

figure5=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure5.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k, label='Hb coefficient')
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel='Hb coeff') 
figure5.savefig('k')

figure5=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure5.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k, label='Hb coefficient')
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel='Hb coeff') 
figure5.savefig('k')

figure5=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure5.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k, label='Hb coefficient')
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel='Hb coeff') 
figure5.savefig('k')

figure5=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure5.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k, label='Hb coefficient')
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel='Hb coeff') 
figure5.savefig('k')


        
        
        
        
        
        
        
        
        



#print("Splitter to Alveoli: %f \nSat Alveolar O2: %f \nSat Alveolar CO2: %f \nSat Alveolar H2O: %f \n" % (n_sat_alveoli, n_sat_alveoli*SAT_PCT_O2, n_sat_alveoli*SAT_PCT_CO2, n_sat_alveoli*SAT_PCT_H2O))
#print("Alveoli to Mixer: %f \nAlveolar O2: %f \nAlveolar CO2: %f \nAlveolar H2O: %f \n" % (n_alv, n_alv*ALV_PCT_O2, n_alv*ALV_PCT_CO2, n_alv*ALV_PCT_H2O))