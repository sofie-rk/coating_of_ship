from adaptive_stepsize_solver import EmbeddedExplicitRungeKutta
from adaptive_stepsize_solver import a, b, c, bhat, order

from values import *

import numpy as np
import matplotlib.pyplot as plt

class MODEL:

    def __call__(self, t, z):
        dzp_dt = 1
        dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(z[1]*L_F - z[0]*L_F)

        return np.array([dzp_dt, dzc_dt])


model = MODEL()




# Run method
fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)
tol = 1.0e-5 # Tolerance in adaptive method
n = 100
Nmax = 10**(n+2)





### FIRST, solve from tau = 0 to an upper limit
t0_0, T_0 = 0, 0.95
z0_0 = np.array([0, 10**(-n)])

ts_0, ys_0 = fehlberg(z0_0, t0_0, T_0, model, Nmax, tol)

## SECOND, solve from tau = limit to tau = 1
t0_1, T_1 = ts_0[-1], 1
z0_1 = ys_0[-1]

ts_1, ys_1 = fehlberg(z0_1, t0_1, T_1, model, Nmax, tol)

ts = np.concatenate((ts_0, ts_1))
ys = np.concatenate((ys_0, ys_1))

plt.plot(ts, ys, '.')
plt.legend(["$z_p$", "$z_c$"])
plt.xlabel("Tau")
plt.grid(True)
plt.show()