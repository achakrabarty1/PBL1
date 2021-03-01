import humidifier as hum
import splitter as sp
from functions import volumetric_to_molar
from functions import getToBlood

# data constants
SAT_PCT_O2 = 0.1967
SAT_PCT_CO2 = 0.003
SAT_PCT_H2O = 0.062
ALV_PCT_O2 = 0.136
ALV_PCT_CO2 = 0.053
ALV_PCT_H2O = 0.062

# calculate outlets and inlets
n_toBlood_O2 = calcToBlood()
n_sat_alveoli = sp.getSplitterAlv()
n_fromBlood_CO2 = n_toBlood_O2
n_alv = n_sat_alveoli

print("Splitter to Alveoli: %f \nSat Alveolar O2: %f \nSat Alveolar CO2: %f \nSat Alveolar H2O: %f \n" % (n_sat_alveoli, n_sat_alveoli*SAT_PCT_O2, n_sat_alveoli*SAT_PCT_CO2, n_sat_alveoli*SAT_PCT_H2O))
print("Alveoli to Mixer: %f \nAlveolar O2: %f \nAlveolar CO2: %f \nAlveolar H2O: %f \n" % (n_alv, n_alv*ALV_PCT_O2, n_alv*ALV_PCT_CO2, n_alv*ALV_PCT_H2O))