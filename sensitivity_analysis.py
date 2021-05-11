#from idealized_route import SENSITIVTY_ANALYSIS_MODEL
from adaptive_stepsize_solver import ODE_solver

from values import L_F, M_unit, rho_p, V_p, rho_c, V_c, M_Cu2O, t_f, M_TBT

from labels import *

from advanced_model import r_TBT_advanced, r_Cu2O_advanced

import numpy as np
import matplotlib.pyplot as plt


class SENSITIVTY_ANALYSIS_MODEL:

    def __init__(self, T, salinity, pH):
        self.T = T
        self.salinity = salinity
        self.pH = pH

    def __call__(self, t, z):

        zp = z[0]
        zc = z[1]

        lc = zc*L_F
        lp = zp*L_F

        r_TBT = r_TBT_advanced(self.T, self.pH, self.salinity)
        r_Cu2O = r_Cu2O_advanced(self.T, self.pH, self.salinity, lc, lp)

        tf = t_f(r_TBT)

        dzpdt = tf * (M_unit * r_TBT) / (L_F * M_TBT * rho_p * V_p)

        dzcdt = tf * (M_Cu2O*r_Cu2O)/(2*V_c*rho_c*L_F) 

        return np.array([dzpdt, dzcdt])

def pH_vary(pH):
    return pH

def sal_vary(sal):
    return sal
    
def temp_vary(T):
    return T 


def varying_temperature():

    temperature = [10+273, 15+273, 20+273, 25+273]

    pH_fixed = 8.2
    sal_fixed = 35.1

    fig, (Xp_ax, Xc_ax) = plt.subplots(1,2)

    for temp in temperature:
        tol = 10e-15
        n = 100
        Nmax = 10**(n+30)

        sensitivity_model = SENSITIVTY_ANALYSIS_MODEL(temp, sal_fixed, pH_fixed)

        t0_s, T_end_s = 0, 0.999

        z0_s = np.array([0, 10**(-n)])

        t_s, y_s = ODE_solver(z0_s, t0_s, T_end_s, sensitivity_model, Nmax, tol)

        r_TBT_c = r_TBT_advanced(temp, pH_fixed, sal_fixed)
        tf_c = t_f(r_TBT_c)

        X_p = np.zeros(len(t_s))
        X_c = np.zeros(len(t_s))
        for i in range(len(t_s)):
            X_p[i] = y_s[i][0]
            X_c[i] = y_s[i][1]

        Xp_ax.plot(t_s*tf_c, X_p, label=temp_label(temp))
        Xc_ax.plot(t_s*tf_c, X_c, label=temp_label(temp))
    
    Xp_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xc_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xp_ax.set_title(Xp_label)
    Xc_ax.set_title(Xc_label)
    Xp_ax.legend(loc="lower right")
    Xc_ax.legend(loc="lower right")
    Xp_ax.grid(True)
    Xc_ax.grid(True)
    fig.suptitle("With constant salinity = " + str(sal_fixed) + " g salt/kg seawater and pH = " + str(pH_fixed))
    plt.show()

def varying_salinity():

    print("STARTING SENSITIVITY ANALYSIS - SALANITY")

    salinities = [35, 40, 45]

    pH_fixed = 8.2
    temp_fixed = 15+273

    fig, (Xp_ax, Xc_ax) = plt.subplots(1,2)

    for sal in salinities:
        print("sal = ", sal)
        tol = 10e-15
        n = 100
        Nmax = 10**(n+10)

        sensitivity_model = SENSITIVTY_ANALYSIS_MODEL(temp_fixed, sal, pH_fixed)

        t0_s, T_end_s = 0, 0.999

        z0_s = np.array([0, 10**(-n)])

        t_s, y_s = ODE_solver(z0_s, t0_s, T_end_s, sensitivity_model, Nmax, tol)

        r_TBT_c = r_TBT_advanced(temp_fixed, pH_fixed, sal)
        tf_c = t_f(r_TBT_c)

        X_p = np.zeros(len(t_s))
        X_c = np.zeros(len(t_s))
        for i in range(len(t_s)):
            X_p[i] = y_s[i][0]
            X_c[i] = y_s[i][1]

        Xp_ax.plot(t_s*tf_c, X_p, label=salinity_label(sal))
        Xc_ax.plot(t_s*tf_c, X_c, label=salinity_label(sal))
    
    Xp_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xc_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xp_ax.set_title(Xp_label)
    Xc_ax.set_title(Xc_label)
    Xp_ax.legend(loc="lower right")
    Xc_ax.legend(loc="lower right")
    Xp_ax.grid(True)
    Xc_ax.grid(True)
    fig.suptitle("With constant T = " + str(temp_fixed-273) + "$^o$C and pH = " + str(pH_fixed))
    plt.show()


def varying_pH():

    print("SENSITIVTY ANALYSIS - pH")

    pHs = [7, 7.5, 8, 8.5]

    temp_fixed = 15+273
    sal_fixed = 35.1

    fig, (Xp_ax, Xc_ax) = plt.subplots(1,2)

    for pH_v in pHs:
        tol = 10e-15
        n = 100
        Nmax = 10**(n+50)

        sensitivity_model = SENSITIVTY_ANALYSIS_MODEL(temp_fixed, sal_fixed, pH_v)

        t0_s, T_end_s = 0, 0.999

        z0_s = np.array([0, 10**(-n)])

        t_s, y_s = ODE_solver(z0_s, t0_s, T_end_s, sensitivity_model, Nmax, tol)

        r_TBT_c = r_TBT_advanced(temp_fixed, pH_v, sal_fixed)
        tf_c = t_f(r_TBT_c)
        print("pH = ", pH_v, "tf_c = ", tf_c)

        X_p = np.zeros(len(t_s))
        X_c = np.zeros(len(t_s))
        for i in range(len(t_s)):
            X_p[i] = y_s[i][0]
            X_c[i] = y_s[i][1]

        Xp_ax.plot(t_s*tf_c, X_p, label=pH_label(pH_v))
        Xc_ax.plot(t_s*tf_c, X_c, label=pH_label(pH_v))
    
    Xp_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xc_ax.set(xlabel=x_label_day, ylabel=conversion_label)
    Xp_ax.set_title(Xp_label)
    Xc_ax.set_title(Xc_label)
    Xp_ax.legend(loc="lower right")
    Xc_ax.legend(loc="lower right")
    Xp_ax.grid(True)
    Xc_ax.grid(True)
    fig.suptitle("With constant T = " + str(temp_fixed-273) + "$^o$C and salinity = " + str(sal_fixed) + " g salt/kg seawater")
    plt.show()

#varying_temperature()
#varying_salinity()
varying_pH()


    