from values import k_L
from values import k1, k2, k_Cu2O, k_1
from values import rho_seawater, Mm_Cl
from values import M_unit, rho_p, V_p
from values import M_Cu2O, V_c, rho_c
from values import D_CuCl

from adaptive_stepsize_solver import EmbeddedExplicitRungeKutta
from adaptive_stepsize_solver import a, b, c, bhat, order

import numpy as np
import matplotlib.pyplot as plt


class ADVANCED_MODEL:

    def __init__(self, temperature, pH, salinity):
        self.temperature = temperature
        self.pH = pH
        self.salinity = salinity

    def __call__(self, t, y):
        '''
        T: temperature [K]
        pH: pH
        salinity: [g salt / kg seawater]
        '''

        lp = y[0]
        lc = y[1]
        
        conc_OH = (10**(-14 + pH))*1000                     # [mol/m3]
        conc_Cl = 0.55 * salinity * rho_seawater / Mm_Cl    # [mol/m3]
        conc_H = 10**(-pH)*1000                             # [mol/m3]
        
        conc_CuCl = 1 / (1/2*k_1 + k_L*D_CuCl/(k_L*(lc-lp) + D_CuCl)) * k_Cu2O(temperature) * conc_H * conc_Cl**2

        r_TBT = (k1(temperature)*conc_OH**0.32*conc_Cl) / (1 + k2*conc_OH**0.43)

        r_Cu2O = k_Cu2O(T)*conc_H*conc_Cl**2 - 1/2 * k_1 * conc_CuCl
        
        dlp_dt = r_TBT * M_unit / (rho_p * V_p)
        dzc_dt = M_Cu2O/(2*V_c*rho_c) * r_Cu2O

        return np.array([dlp_dt, dzc_dt])






# Run method
fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)
tol = 10e-5
n = 100
Nmax = 10**(n+3)

### SOLVE FROM t = 0 to t = 400
temperature = 10 + 273                    # [K]
pH = 8.2
salinity = 35.1                 # [g salt / kg seawater]

model = ADVANCED_MODEL(temperature, pH, salinity)

t0, T = 0, 400

l0 = np.array([0, 10**(-n)])

ts, ys = fehlberg(l0, t0, T, model, Nmax, tol)

thickness = np.zeros(len(ts))
for i in range(len(ts)):
    thickness[i] = ys[i][1] - ys[i][0]

plt.plot(ts, ys, '.')
plt.legend(["$l_p$", "$l_c$"])
plt.xlabel("Time t [days]")
plt.ylabel("Length [m]")
plt.show()

plt.plot(ts, thickness)
plt.xlabel("Time t[days]")
plt.ylabel("Thickness [m]")
plt.title("Thickness of the leached layer")
plt.show()
