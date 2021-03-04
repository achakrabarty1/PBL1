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
nHmc = np.zeros(l)
nHw = np.zeros(l)
nHOc = np.zeros(l)

ngc = np.zeros(l)
nCg = np.zeros(l)
nHfg = np.zeros(l)
nH2Og = np.zeros(l)
nH20mc = np.zeros(l)
ngmc = np.zeros(l)
nHfmc = np.zeros(l)
'''
nHE[0] = n_Hb_EPO
nHcm[0] = n_Hb_cap_met
nHOcm[0] = n_HbO_cap_met
nHfcm[0] = n_Hbf_cap_met    
'''   
    
for i in range(-1,len(a)-1):
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
    
    nHE[i] = (k[i]-1)*normal_n_Hb_EPO
    nHcm[i] = Dens_Hb*CO/MW_Hb+nHE[i]

    nHOcm[i] = nOac[i]/4
    nHfcm[i] = nHcm[i]-nHOcm[i]+nHE[i]
    nC = .0225*V_b   #based on volumetric flow rate and conc
    nHmc[i]=nHfcm[i]+nHOcm[i]
    nHw[i] = nHmc[i] + nHE[i] - nHcm[i]
    ngi = V_b * c_gluc_in / mw_gluc
    nH2Oi = m_H2O_in * 24 * 60 / mw_H2O
    ngc[i]=nOc[i]/6
    nHOc[i] = nOc[i] / 4
    

# calculations of generated products
# because we found O2 to be the limiting reagent, and RQ = 1:
    nCg[i] = nOc[i]
    nH2Og[i] = nOc[i]
    nHfg[i] = nHOc[i]
    nH20mc[i] = nH2Oi+nH2Og[i]
    ngmc[i] = ngi-ngc[i]
    nHfmc[i] = nHfcm[i] + nHfg[i]

 
figure=plt.figure(num=1, clear=True)  #plotting for healthy v arteriosclerosis
axis=figure.add_subplot(1,1,1)
axis.grid(True)    
axis.plot(a, vtbO)
axis.legend()        
axis.set(xlabel='Area (cm^2)', ylabel = 'O2 flow (L/min)') 
title1=plt.title('O2 flow from alveoli to capillaries')
figure.savefig('vtbo')  
        
figure2=plt.figure(num=2, clear=True) #plotting for arteriosclerosis v treated
axis2=figure2.add_subplot(1,1,1) 
axis2.grid(True)     
axis2.plot(a, ntbO) 
title2=plt.title('Molar O2 flow from Alv to Cap')      
axis2.legend()
axis2.set(xlabel='Area (cm^2)', ylabel='O2 flow (mol/min)') 
figure2.savefig('ntbo')  

figure3=plt.figure(num=3, clear=True)  #plotting for healthy v arteriosclerosis
axis3=figure3.add_subplot(1,1,1)
axis3.grid(True)    
axis3.plot(a, vtbC)
axis3.legend()        
axis3.set(xlabel='Area (cm^2)', ylabel = 'CO2 flow (L/min)') 
title3=plt.title('CO2 flow from Alv to Cap')
figure3.savefig('vtbc')

figure4=plt.figure(num=4, clear=True)  #plotting for healthy v arteriosclerosis
axis4=figure4.add_subplot(1,1,1)
axis4.grid(True)    
axis4.plot(a, ntbC)
axis4.legend()        
axis4.set(xlabel='Area (cm^2)', ylabel = 'CO2 flow (mol/min)') 
title4=plt.title('CO2 flow from alveoli to capillaries')
figure4.savefig('ntbc')

figure5=plt.figure(num=5, clear=True)  #plotting for healthy v arteriosclerosis
axis5=figure5.add_subplot(1,1,1)
axis5.grid(True)    
axis5.plot(a, k)
axis5.legend()        
axis5.set(xlabel='Area (cm^2)', ylabel = 'Hemoglobin coefficient (dimensionless)') 
title5=plt.title('Hemoglobin coefficient')
figure5.savefig('k')

figure6=plt.figure(num=6, clear=True)  #plotting for healthy v arteriosclerosis
axis6=figure6.add_subplot(1,1,1)
axis6.grid(True)    
axis6.plot(a, nOmc)
title6=plt.title('Oxygen flow from Metabolic Reactor to Capillaries')
axis6.legend()        
axis6.set(xlabel='Area (cm^2)', ylabel='O2 (mol/min)') 
figure6.savefig('nOmc')

figure7=plt.figure(num=7, clear=True)  #plotting for healthy v arteriosclerosis
axis7=figure7.add_subplot(1,1,1)
axis7.grid(True)    
axis7.plot(a, k)
title7=plt.title('O2 flow from capillaries to metabolic reactor')
axis7.legend()        
axis7.set(xlabel='Area (cm^2)', ylabel='O2 flow (mol/min)') 
figure7.savefig('nOcm')

figure8=plt.figure(num=8, clear=True)  #plotting for healthy v arteriosclerosis
axis8=figure8.add_subplot(1,1,1)
axis8.grid(True)    
axis8.plot(a, nOc)
axis8.legend()        
axis8.set(xlabel='Area (cm^2)', ylabel='Molar flow (mol/min)') 
title8=plt.title('Oxygen consumed in metabolic reactor')
figure8.savefig('nOc')

figure9=plt.figure(num=9, clear=True)  #plotting for healthy v arteriosclerosis
axis9=figure9.add_subplot(1,1,1)
axis9.grid(True)    
axis9.plot(a, nCac)
axis9.legend() 
title9=plt.title('CO2 diffusion from alveoli to capillaries')       
axis9.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure9.savefig('nCac')

figure10=plt.figure(num=10, clear=True)  #plotting for healthy v arteriosclerosis
axis10=figure10.add_subplot(1,1,1)
axis10.grid(True)    
axis10.plot(a, k)
axis10.legend()        
axis10.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
title10=plt.title('CO2 diffusion from capillaries to alveoli')
figure10.savefig('nCca')

figure11=plt.figure(num=11, clear=True)  #plotting for healthy v arteriosclerosis
axis11=figure11.add_subplot(1,1,1)
axis11.grid(True)    
axis11.plot(a, nHE)
axis11.legend()     
axis11.set(xlabel='Area (cm^2)', ylabel = 'Molar flow rate (mol/min)') 
title11=plt.title('Hemoglobin contribution from EPO facilitation')
figure11.savefig('nHE')

figure12=plt.figure(num=12, clear=True)  #plotting for healthy v arteriosclerosis
axis12=figure12.add_subplot(1,1,1)
axis12.grid(True)    
axis12.plot(a, nHcm)
axis12.legend()   
title12=plt.title('Hemoglobin flow from capillaries to metabolic reactor')     
axis12.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure12.savefig('nHcm')

figure13=plt.figure(num=13, clear=True)  #plotting for healthy v arteriosclerosis
axis13=figure13.add_subplot(1,1,1)
axis13.grid(True)    
axis13.plot(a, nHOcm)
axis13.legend()   
title13=plt.title('Oxyhemoglobin flow from capillaries to metabolic reactor')     
axis13.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure13.savefig('nHOcm')

figure14=plt.figure(num=14, clear=True)  #plotting for healthy v arteriosclerosis
axis14=figure14.add_subplot(1,1,1)
axis14.grid(True)    
axis14.plot(a, nHfcm)
axis14.legend()   
title14=plt.title('Free Hemoglobin flow from capillaries to metabolic reactor')     
axis14.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure14.savefig('nHfcm')

figure15=plt.figure(num=15, clear=True)  #plotting for healthy v arteriosclerosis
axis15=figure15.add_subplot(1,1,1)
axis15.grid(True)    
axis15.plot(a, nHw)
axis15.legend()   
title15=plt.title('Total hemoglobin flow to waste')     
axis15.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure12.savefig('nHw')

figure16=plt.figure(num=16, clear=True)  #plotting for healthy v arteriosclerosis
axis16=figure16.add_subplot(1,1,1)
axis16.grid(True)    
axis16.plot(a, ngc)
axis16.legend()   
title16=plt.title('Total glucose flow from storage and food consumption')     
axis16.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure16.savefig('ngc')

figure17=plt.figure(num=17, clear=True)  #plotting for healthy v arteriosclerosis
axis17=figure17.add_subplot(1,1,1)
axis17.grid(True)    
axis17.plot(a, nHOc)
axis17.legend()   
title17=plt.title('Oxyhemoglobin consumed for glycolysis')     
axis17.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure17.savefig('nHOc')

figure18=plt.figure(num=18, clear=True)  #plotting for healthy v arteriosclerosis
axis18=figure18.add_subplot(1,1,1)
axis18.grid(True)    
axis18.plot(a, nCg)
axis18.legend()   
title18=plt.title('CO2 generated from glycolysis')     
axis18.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure18.savefig('nCg')

figure19=plt.figure(num=19, clear=True)  #plotting for healthy v arteriosclerosis
axis19=figure19.add_subplot(1,1,1)
axis19.grid(True)    
axis19.plot(a, nHfg)
axis19.legend()   
title19=plt.title('Free hemoglobin generated from glycolysis')     
axis19.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure19.savefig('nHfg')

figure20=plt.figure(num=20, clear=True)  #plotting for healthy v arteriosclerosis
axis20=figure20.add_subplot(1,1,1)
axis20.grid(True)    
axis20.plot(a, nH2Og)
axis20.legend()   
title20=plt.title('H2O) generated from glycolysis')     
axis20.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure20.savefig('nH2Og')

figure21=plt.figure(num=21, clear=True)  #plotting for healthy v arteriosclerosis
axis21=figure21.add_subplot(1,1,1)
axis21.grid(True)    
axis21.plot(a, ngmc)
axis21.legend()   
title21=plt.title('Glucose flow from metabolic reactor to capillaries')     
axis21.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure21.savefig('ngmc')

figure22=plt.figure(num=22, clear=True)  #plotting for healthy v arteriosclerosis
axis22=figure22.add_subplot(1,1,1)
axis22.grid(True)    
axis22.plot(a, nHfmc)
axis22.legend()   
title22=plt.title('Free hemoglobin travelling from metabolic reactor to capillaries')     
axis22.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure22.savefig('nHfmc')

figure23=plt.figure(num=23, clear=True)  #plotting for healthy v arteriosclerosis
axis23=figure23.add_subplot(1,1,1)
axis23.grid(True)    
axis23.plot(a, nCg)
axis23.legend()   
title23=plt.title('CO2 generated from glycolysis')     
axis23.set(xlabel='Area (cm^2)', ylabel='Molar flow rate (mol/min)') 
figure23.savefig('nCg')



        
        
        
        
        
        
        
        
        



#print("Splitter to Alveoli: %f \nSat Alveolar O2: %f \nSat Alveolar CO2: %f \nSat Alveolar H2O: %f \n" % (n_sat_alveoli, n_sat_alveoli*SAT_PCT_O2, n_sat_alveoli*SAT_PCT_CO2, n_sat_alveoli*SAT_PCT_H2O))
#print("Alveoli to Mixer: %f \nAlveolar O2: %f \nAlveolar CO2: %f \nAlveolar H2O: %f \n" % (n_alv, n_alv*ALV_PCT_O2, n_alv*ALV_PCT_CO2, n_alv*ALV_PCT_H2O))