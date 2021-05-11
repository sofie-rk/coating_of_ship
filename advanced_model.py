from values import k_L
from values import k1, k2, k_Cu2O, k_1
from values import rho_seawater, Mm_Cl
from values import M_unit, rho_p, V_p, M_TBT
from values import M_Cu2O, V_c, rho_c
from values import D_CuCl, L_F

from adaptive_stepsize_solver import ODE_solver

from conditions_advanced import *

import numpy as np
import matplotlib.pyplot as plt


def r_TBT_advanced(T, pH, sal):

    conc_OH = OH_conc(pH)
    conc_Cl = Cl_conc(sal)

    r_TBT_a = (k1(T)*conc_OH**0.32*conc_Cl) / (1 + k2*conc_OH**0.43) * M_TBT

    return r_TBT_a

def r_Cu2O_advanced(T, pH, sal, lc, lp):

    conc_H = H_conc(pH)
    conc_Cl = Cl_conc(sal)

    conc_CuCl = 1 / (1/2*k_1 + k_L*D_CuCl/(k_L*(lc-lp) + D_CuCl)) * k_Cu2O(T) * conc_H * conc_Cl**2
    r_Cu2O = k_Cu2O(T)*conc_H*conc_Cl**2 - 1/2 * k_1 * conc_CuCl

    return r_Cu2O




class ADVANCED_MODEL:

    def __init__(self, T, salinity, pH):
        self.T = T
        self.salinity = salinity
        self.pH = pH

    def __call__(self, t, l):

        lp = l[0]
        lc = l[1]

        r_TBT = r_TBT_advanced(self.T(t), self.pH(t), self.salinity(t))
        r_Cu2O = r_Cu2O_advanced(self.T(t), self.pH(t), self.salinity(t), lc, lp)

        dlpdt = (M_unit * r_TBT) / (M_TBT * rho_p * V_p)

        dlcdt = (M_Cu2O*r_Cu2O)/(2*V_c*rho_c) 

        return np.array([dlpdt, dlcdt])





# Run method
tol = 10e-20
n = 100
Nmax = 10**(n+10)

### FIRST 400 days, Bueno Aires coastline ### 

model_400 = ADVANCED_MODEL(p_400_temp, p_400_sal, p_400_pH)

t0_400, T_end_400 = 0, 400

l0_400 = np.array([0, 10**(-n)])

t_400, y_400 = ODE_solver(l0_400, t0_400, T_end_400, model_400, Nmax, tol)


### 20 DAYS VOYAGE ###
model_20 = ADVANCED_MODEL(p_20_temp, p_20_sal, p_20_pH)

t0_20, T_end_20 = 0, 19.9

l0_20 = np.array([y_400[-1][0], y_400[-1][1]])

t_20, y_20 = ODE_solver(l0_20, t0_20, T_end_20, model_20, Nmax, tol)


### 400 + 20 days modeling
t_1 = t_400[:]
y_1 = y_400[:]

t_2 = t_20[:] + 400
y_2 = y_20[:]

t_420 = np.concatenate((t_1, t_2))
y_420 = np.concatenate((y_1, y_2))




