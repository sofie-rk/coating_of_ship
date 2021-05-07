from values import k_L
from values import k1, k2, k_Cu2O, k_1
from values import rho_seawater, Mm_Cl
from values import M_unit, rho_p, V_p, M_TBT
from values import M_Cu2O, V_c, rho_c
from values import D_CuCl, L_F

from adaptive_stepsize_solver import EmbeddedExplicitRungeKutta
from adaptive_stepsize_solver import a, b, c, bhat, order

import numpy as np
import matplotlib.pyplot as plt


class ADVANCED_MODEL:

    def __call__(self, t, l):

        lp = l[0]
        lc = l[1]

        conc_OH = (10**(-14 + pH))*1000                     # [mol/m3]
        conc_Cl = 0.55 * salinity * rho_seawater / Mm_Cl    # [mol/m3]
        conc_H = 10**(-pH)*1000                             # [mol/m3]

        k1_i = k1(T)
        k2_i = k2


        r_TBT = (k1(T)*conc_OH**0.32*conc_Cl) / (1 + k2*conc_OH**0.43) * M_TBT


        
        conc_CuCl = 1 / (1/2*k_1 + k_L*D_CuCl/(k_L*(lc-lp) + D_CuCl)) * k_Cu2O(T) * conc_H * conc_Cl**2

        r_Cu2O = k_Cu2O(T)*conc_H*conc_Cl**2 - 1/2 * k_1 * conc_CuCl

        dlpdt = (M_unit * r_TBT) / (M_TBT * rho_p * V_p)

        dlcdt = (M_Cu2O*r_Cu2O)/(2*V_c*rho_c) 

        return np.array([dlpdt, dlcdt])

        
# Run method
fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)
tol = 10e-12
n = 100
Nmax = 10**(n+7)

### SOLVE FROM t = 0 to t = 400
T = 15.94 + 273                    # [K]
pH = 8.2
salinity = 35.46                 # [g salt / kg seawater]

model = ADVANCED_MODEL()

t0, T_end = 0, 400

z0 = np.array([0, 10**(-n)])

ts, ys = fehlberg(z0, t0, T_end, model, Nmax, tol)

print("z_p, 400 days: ", ys[-1][0])
print("z_c, 400 days: ", ys[-1][1])

thickness = np.zeros(len(ts))
for i in range(len(ts)):
    thickness[i] = ys[i][1] - ys[i][0]

plt.plot(ts, ys*10**3)
plt.legend(["$l_p$", "$l_c$"])
plt.grid(True)
plt.ylabel("[mm]")
plt.xlabel("Time [days]")
plt.show()

plt.plot(ts, thickness*10**3)
plt.xlabel("Time t[days]")
plt.ylabel("Thickness [mm]")
plt.title("Thickness of the leached layer")
plt.grid(True)
plt.show()
