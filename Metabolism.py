#import alveoli
import alveoli as alv
import capillaries as cap

# basis for calculations
mw_gluc = 180.16 #g/mol
c_gluc_in = 0.9 #g/L
V_b = 5.0 #L/min
m_H2O_in = 2100 #mL/day
mw_H2O = 18.015 #g/mol

# calculations of intake
n_O2_in = alv.getO2() #mol/min
n_HbO_in = n_O2_in / 4
n_Hb_in = cap.getHb()
n_Hbf_in = n_Hb_in - n_HbO_in
n_gluc_in = V_b * c_gluc_in / mw_gluc #mol/min
n_H2O_in = m_H2O_in * 24 * 60 / mw_H2O

# determination of limiting reagent
LR = min(n_O2_in/6, n_gluc_in)
print(LR)

# calculations of consumed products
# because we found O2 to be the limiting reagent,
n_O2_con = n_O2_in
n_gluc_con = n_O2_con / 6
n_HbO_con = n_O2_con / 4

# calculations of generated products
# because we found O2 to be the limiting reagent, and RQ = 1:
n_CO2_gen = n_O2_in
n_H2O_gen = n_O2_in
n_Hbf_gen = n_HbO_con

# calculations of output
n_CO2_out = n_CO2_gen
n_H2O_out = n_H2O_in + n_H2O_gen
n_O2_out = n_O2_in - n_O2_con
n_gluc_out = n_gluc_in - n_gluc_con
n_Hbf_out = n_Hbf_in + n_Hbf_gen
n_HbO_out = n_HbO_in - n_HbO_con

def get_CO2_met():
  return n_CO2_out
  
def get_H2O_met():
  return n_H2O_out

def get_Hbf_met():
  return n_Hbf_out