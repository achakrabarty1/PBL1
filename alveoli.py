# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 01:30:54 2021

@author: stper
"""

import humidifier as hum
import splitter as sp
from functions import volumetric_to_molar
from functions import getToBlood

def getGradient(pBlood, pAlv):
    return (pBlood-pAlv)/DIST
def calcToBlood(coeff,gradient):
    return (coeff)*(ALV_AREA)*(gradient)
def getNBloodO2():
    return n_toBlood_O2
def getVBloodO2():
    return v_toBlood_O2
def getNBloodC2():
    return n_toBlood_CO2
def getVBloodCO2():
    return v_toBlood_CO2

# data constants
SAT_PCT_O2 = 0.1967
SAT_PCT_CO2 = 0.003
SAT_PCT_H2O = 0.062
ALV_PCT_O2 = 0.136
ALV_PCT_CO2 = 0.053
ALV_PCT_H2O = 0.062

ALV_AREA = 1.4 * (10^6)  #cm^2
PO2_BLOOD = 104 #mm hg
PO2_ALV = 40 #mm hg
PO2_GRADIENT = getGradient(PO2_BLOOD, PO2_ALV)
DIFF_COEFF_O2 = 2.4*(10^-5) #cm^2/s

PCO2_BLOOD = 40
PCO2_ALV = 45
PCO2_GRADIENT = getGradient(PCO2_BLOOD, PCO2_ALV)
DIFF_COEFF_CO2 = 1.55*(10^-5) #cm^2/s

DIST = .03 #distance from alveola to capillary


# calculate outlets and inlets
v_toBlood_O2 = calcToBlood(DIFF_COEFF_O2, PO2_GRADIENT)
n_toBlood_O2 = volumetric_to_molar(v_toBlood_O2)
v_toBlood_CO2 = calcToBlood(DIFF_COEFF_CO2, PCO2_GRADIENT)
n_toBlood_CO2 = volumetric_to_molar(v_toBlood_CO2)
n_sat_alveoli = sp.getSplitterAlv()
n_fromBlood_CO2 = n_toBlood_O2
n_alv = n_sat_alveoli





print("Splitter to Alveoli: %f \nSat Alveolar O2: %f \nSat Alveolar CO2: %f \nSat Alveolar H2O: %f \n" % (n_sat_alveoli, n_sat_alveoli*SAT_PCT_O2, n_sat_alveoli*SAT_PCT_CO2, n_sat_alveoli*SAT_PCT_H2O))
print("Alveoli to Mixer: %f \nAlveolar O2: %f \nAlveolar CO2: %f \nAlveolar H2O: %f \n" % (n_alv, n_alv*ALV_PCT_O2, n_alv*ALV_PCT_CO2, n_alv*ALV_PCT_H2O))