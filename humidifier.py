from functions import volumetric_to_molar

# data values basis for calculation
TIDAL_VOLUME = 6000 # ml/min
AIR_PCT_O2 = 0.21
AIR_PCT_CO2 = 0.0004
AIR_PCT_H2O = 0.005
SAT_PCT_H2O = 0.062

# calculate molar flow rates of air stream
n_air_tot = volumetric_to_molar(TIDAL_VOLUME)
n_air_O2 = n_air_tot * AIR_PCT_O2
n_air_CO2 = n_air_tot * AIR_PCT_CO2
n_air_H2O = n_air_tot * AIR_PCT_H2O

# calculate molar flow rate in URT stream using equation derived
n_URT_H2O = n_air_tot * (SAT_PCT_H2O - AIR_PCT_H2O) / (1 - SAT_PCT_H2O)

# calculate molar flow rates in saturated stream
# O2 and CO2 unchanged from air stream
n_sat_O2 = n_air_O2
n_sat_CO2 = n_air_CO2
n_sat_H2O = n_air_H2O + n_URT_H2O
n_sat_tot = n_air_tot + n_URT_H2O


print("Air Stream (mol/min) \nTotal: %f \nOxygen: %f \nCO2: %f \nWater: %f \n" % (n_air_tot, n_air_O2, n_air_CO2, n_air_H2O))
print("URT Water (mol/min): %f\n" % n_URT_H2O)
print("Saturated Stream (mol/min) \nTotal: %f \nOxygen: %f \nCO2: %f \nWater: %f \n" % (n_sat_tot, n_sat_O2, n_sat_CO2, n_sat_H2O))

def getTot():
    return n_sat_tot
def getO2():
    return n_sat_O2
def getCO2():
    return n_sat_CO2
def getH2O():
    return n_sat_H2O