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
tol = 1.0e-15 # Tolerance in adaptive method
n = 100
Nmax = 10**(n+8)





### FIRST, solve from tau = 0 to an upper limit
t0_0, T_0 = 0, 0.8
z0_0 = np.array([0, 10**(-n)])

ts_0, ys_0 = fehlberg(z0_0, t0_0, T_0, model, Nmax, tol)



## SECOND, solve from tau = limit to tau = 1
t0_1, T_1 = ts_0[-1], 1
z0_1 = ys_0[-1]

ts_1, ys_1 = fehlberg(z0_1, t0_1, T_1, model, Nmax, tol)

tau_simple = np.concatenate((ts_0, ts_1))
X_simple = np.concatenate((ys_0, ys_1))

def plot_conversion_simple():
    plt.plot(tau_simple*tf, X_simple)
    plt.legend(["$X_p$", "$X_c$"])
    plt.ylabel("Conversion [-]")
    plt.xlabel("Time [days]")
    plt.title("Simple model, conversion vs time. Dimensionless model used.")
    plt.grid(True)
    plt.show()


### THICKNESS OF THE LEACHED LAYER ###

def plot_thickness_simple():
    thickness = np.zeros(len(tau_simple))
    for i in range(len(tau_simple)):
        thickness[i] = (X_simple[i][1] - X_simple[i][0]) * L_F * 10**3  #[mm]
    plt.plot(tau_simple[:-1]*tf, thickness[:-1])
    plt.title("Thickness of leached layer. Simple model. Dimensionless model used.")
    plt.xlabel("Time [days]")
    plt.ylabel("$l_c - l_p$ [mm]")
    plt.show()


### RELEASE OF Cu2+

def plot_release_Cu():

    r_Cu2O_simple = np.zeros(len(tau_simple))

    for i in range(len(r_Cu2O_simple)):

        r_Cu2O_simple[i] = D_CuCl*C_CuCl_s/(L_F*(X_simple[i][1] - X_simple[i][0])) * M_Cu2O * 10**(9) * 10**(-4) # [micro g/cm2 day]
    
    plt.plot(tau_simple[1500:-10]*tf, r_Cu2O_simple[1500:-10])
    plt.xlabel("Time [days]")
    plt.ylabel(r"Release of $Cu^{2+}$  [$\mu$ g / $cm^2$ day]")
    plt.show()



plot_release_Cu()