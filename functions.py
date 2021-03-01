# STP: Pressure = 760 mmHg, Temp = 273 K (0 C), R = 62.363 L mmHg mol-1 K-1 = 6.2363e4 cm^3 mmHg mol-1 K-1
STP_P = 760
STP_T = 273
R = 6.2363e4
# Converts volumetric flow rate to molar flow rate using Ideal Gas Law
# Arguments: temperature, pressure, volumetric flow rate in mL/min
# from Ideal Gas Law PV = nRT -> n = (P/RT)V, calculate how many moles of air are in 1000mL gas
# default is STP
def volumetric_to_molar(vol_flow, pressure = STP_P, temperature = STP_T):
    conv_factor = pressure / (R * temperature)
    return vol_flow * conv_factor
