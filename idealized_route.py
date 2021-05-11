from values import *
from conditions_advanced import *
from labels import *

from advanced_model import r_TBT_advanced, r_Cu2O_advanced

from adaptive_stepsize_solver import ODE_solver

import numpy as np
import matplotlib.pyplot as plt


# Choosing tf average to be equal to the simple tf
tf_average = tf_simple


class IDEALIZED_ROUTE:

    def __init__(self, T, salinity, pH):
        self.T = T
        self.salinity = salinity
        self.pH = pH

    def __call__(self, t, z):

        zp = z[0]
        zc = z[1]

        lc = zc*L_F
        lp = zp*L_F

        r_TBT = r_TBT_advanced(self.T(t), self.pH(t), self.salinity(t))
        r_Cu2O = r_Cu2O_advanced(self.T(t), self.pH(t), self.salinity(t), lc, lp)

        dzpdt = tf_average * (M_unit * r_TBT) / (L_F * M_TBT * rho_p * V_p)

        dzcdt = tf_average * (M_Cu2O*r_Cu2O)/(2*V_c*rho_c*L_F) 

        return np.array([dzpdt, dzcdt])

# Run method
tol = 10e-10
n = 100
Nmax = 10**(n+15)

idealized_route_model = IDEALIZED_ROUTE(p_i_temp, p_i_sal, p_i_pH)

t0_s, T_end_s = 0, 0.999

z0_s = np.array([0, 10**(-n)])

t_s, y_s = ODE_solver(z0_s, t0_s, T_end_s, idealized_route_model, Nmax, tol)

def plot_conversion_idealized():

    plt.plot(t_s, y_s)
    plt.plot([0.375, 0.375],[0, 0.7], '--', color="black")
    plt.plot([0.500, 0.500],[0, 0.7], '--', color="black")
    plt.plot([0.625, 0.625],[0, 0.7], '--', color="black")
    plt.plot([0.750, 0.750],[0, 0.7], '--', color="black")
    plt.plot([0.875, 0.875],[0, 0.7], '--', color="black")
    plt.legend(conversion_legend)
    plt.grid(True)
    plt.ylim(0, 0.75)
    plt.xlabel(tau_label)
    plt.ylabel(conversion_label)
    plt.show()

def plot_release_Cu():
    r_Cu = np.zeros(len(t_s))

    for i in range(len(r_Cu)):

        t = t_s[i]

        r_Cu[i] = r_Cu2O_advanced(p_i_temp(t), p_i_pH(t), p_i_sal(t), y_s[i][1]*L_F, y_s[i][0]*L_F) * M_Cu * 10**(-4)

    plt.plot(t_s, r_Cu)
    plt.plot([0.375, 0.375],[0, 30], '--', color="black")
    plt.plot([0.500, 0.500],[0, 30], '--', color="black")
    plt.plot([0.625, 0.625],[0, 30], '--', color="black")
    plt.plot([0.750, 0.750],[0, 30], '--', color="black")
    plt.plot([0.875, 0.875],[0, 30], '--', color="black")
    plt.xlabel("Tau [-]")
    plt.ylabel(release_Cu_label)
    plt.grid(True)
    plt.ylim(0,40)
    plt.show()

plot_conversion_idealized()
plot_release_Cu()
