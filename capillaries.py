# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 01:19:39 2021

@author: stper
"""

import humidifier as hum
import splitter as sp
import metabolism as met
import alveoli as alv
from functions import volumetric_to_molar
from functions import getToBlood

pp_O2 = 0.197     #not likely to be used, as we care about bound O2
pp_CO2 = .09

k = alv.getNormalArea()/alv.getArea()   #add another constant here? Not sure how this would scale
                                        #as i'm not sure how much SA is lost with emphysema

MW_Hb = 64500
CO = 5   #5 L blood/min
Dens_Hb = 150 #g/L blood
HbOrel = 0.59  #ideal L O2/(100 g Hb * min)

v_O2_in=alv.getVBloodO2()  #L/min
n_O2_in = alv.getNBloodO2()*k
n_O2_out = n_O2_in  #mol/min

n_Hb_ideal = (1/(HbOrel*100)/MW_Hb)*v_O2_in
n_Hb = Dens_Hb*CO/MW_Hb
n_HbO = n_O2_out/4
n_HbF = n_Hb - n_HbO

n_CO2 = .0225*5   #mol/min
n_CO2_in = alv.getNBloodCO2() 
n_CO2_met = n_O2_out
n_CO2_out = n_CO2_in + n_CO2_met



n_Hb_met_out = n_Hb
n_Hb_met_in = n_Hb
normal_n_Hb_EPO = Dens_Hb*CO/(100*MW_Hb)
n_Hb_EPO = (k)*normal_n_Hb_EPO
n_Hb_waste = n_Hb_met_in + n_Hb_EPO - n_Hb_met_out



def getHb():
    return n_Hb
def getO2():
    return n_O2_out
