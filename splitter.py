from functions import getToBlood
import humidifier as hum
# splitter code
# get alveoli to blood O2 flow rate
SAT_PCT_O2 = 0.1967
SAT_PCT_CO2 = 0.003
SAT_PCT_H2O = 0.062
ALV_PCT_O2 = 0.136

n_toBlood_O2 = calcToBlood()
n_sat_tot = hum.getTot()

n_sat_alveoli = n_toBlood_O2 / (SAT_PCT_O2 - ALV_PCT_O2)
n_sat_deadspace = n_sat_tot - n_sat_alveoli

print("Splitter to Alveoli: %f \nSat Alveolar O2: %f \nSat Alveolar CO2: %f \nSat Alveolar H2O: %f \n" % (n_sat_alveoli, n_sat_alveoli*SAT_PCT_O2, n_sat_alveoli*SAT_PCT_CO2, n_sat_alveoli*SAT_PCT_H2O))
print("Splitter to Deadspace: %f \nDS O2: %f \nDS CO2: %f \nDS H2O: %f \n" % (n_sat_deadspace, n_sat_deadspace*SAT_PCT_O2, n_sat_deadspace*SAT_PCT_CO2, n_sat_deadspace*SAT_PCT_H2O))

def getSplitterAlv():
    return n_sat_alveoli