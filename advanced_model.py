from values import k_L
from values import k1, k2, k_Cu2O, k_1

def advanced_model(t, l, T, pH, salinity):

    conc_OH = 0
    conc_Cl = 0
    conc_H = 0

    conc_CuCl = 0

    r_TBT = (k1(T)*conc_OH**0.32*conc_Cl) / (1 + k2*conc_OH**0.43)

    r_Cu2O = k_Cu2O(T)*conc_H*conc_Cl - 1/2 * k_1 * conc_CuCl