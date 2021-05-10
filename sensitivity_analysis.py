from values import *

from advanced_model import r_TBT_advanced, r_Cu2O_advanced

from adaptive_stepsize_solver import EmbeddedExplicitRungeKutta
from adaptive_stepsize_solver import a, b, c, bhat, order

import numpy as np
import matplotlib.pyplot as plt


def T_s(tau):
    if (tau<0.375):
        return 10+273
    elif (tau>=0.375 and tau<0.5):
        return 15+273
    elif (tau>=0.5 and tau<=1):
        return 20+273

def sal_s(tau):
    if (tau<0.625):
        return 35.1
    elif (tau>=0.625 and tau<0.750):
        return 40.1
    elif (tau>=0.750 and tau<=1):
        return 45.1


def pH_s(tau):
    if (tau<0.875):
        return 8.2
    elif (tau>=0.875 and tau<1):
        return 7.8
    elif (tau == 1):
        return 7.4


def t_f(r_TBT_i):
    #print("r_TBT_i: ", r_TBT_i)
    return (L_F*M_TBT*rho_p*V_p)/(r_TBT_i * M_unit)


def find_tf_average():

    tau = [0, 0.125, 0.250, 0.375, 0.5, 0.625, 0.75, 0.875, 1]
    temp = [10+273, 10+273, 10+273, 15+273, 20+273, 20+273, 20+273, 20+273, 20+273]
    sal = [35.1, 35.1, 35.1, 35.1, 35.1, 40.1, 45.1, 45.1, 45.1]
    pH = [8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 7.8, 7.4]

    tf = 0
    for i in range(len(tau)):
        r_TBT_i = r_TBT_advanced(temp[i], pH[i], sal[i])
        print("t_f: ", t_f(r_TBT_i))
        tf += t_f(r_TBT_i)
    tf = tf/len(tau)
    print("Average t_F: ", tf)

    return tf


tf_average = find_tf_average() - 5000

class SENSITIVTY_ANALYSIS_MODEL:

    def __init__(self, T, salinity, pH):
        self.T = T
        self.salinity = salinity
        self.pH = pH

    def __call__(self, t, l):

        zp = l[0]
        zc = l[1]

        lc = zc*L_F
        lp = zp*L_F

        r_TBT = r_TBT_advanced(self.T(t), self.pH(t), self.salinity(t))
        r_Cu2O = r_Cu2O_advanced(self.T(t), self.pH(t), self.salinity(t), lc, lp)

        dzpdt = tf_average * (M_unit * r_TBT) / (L_F * M_TBT * rho_p * V_p)

        dzcdt = tf_average * (M_Cu2O*r_Cu2O)/(2*V_c*rho_c) 

        return np.array([dzpdt, dzcdt])

# Run method
fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)
tol = 10e-10
n = 100
Nmax = 10**(n+10)

sensitivity_model = SENSITIVTY_ANALYSIS_MODEL(T_s, sal_s, pH_s)

t0_s, T_end_s = 0, 0.999

z0_s = np.array([0, 10**(-n)])

t_s, y_s = fehlberg(z0_s, t0_s, T_end_s, sensitivity_model, Nmax, tol)

plt.plot(t_s, y_s)
# plt.plot([0.375, 0.375],[0, 1], color="black")
# plt.plot([0.500, 0.500],[0, 1], color="black")
# plt.plot([0.625, 0.625],[0, 1], color="black")
# plt.plot([0.750, 0.750],[0, 1], color="black")
# plt.plot([0.875, 0.875],[0, 1], color="black")
plt.legend(["$X_p$", "$X_c$"])
plt.grid(True)
plt.xlabel("Tau [-]")
plt.ylabel("Conversion [-]")
plt.show()