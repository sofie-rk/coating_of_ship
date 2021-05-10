from values import k_L
from values import k1, k2, k_Cu2O, k_1
from values import rho_seawater, Mm_Cl
from values import M_unit, rho_p, V_p, M_TBT
from values import M_Cu2O, V_c, rho_c
from values import D_CuCl, L_F

from adaptive_stepsize_solver import EmbeddedExplicitRungeKutta
from adaptive_stepsize_solver import a, b, c, bhat, order

from interpolation import p_temperature, p_salinity, p_pH
from interpolation import temperatures, salinity, pH

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

def OH_conc(pH):
    return (10**(-14 + pH))*1000    # [mol/m3]

def H_conc(pH):
    return 10**(-pH)*1000           # [mol/m3]

def Cl_conc(sal):
    return 0.55 * sal * rho_seawater / Mm_Cl    # [mol/m3]


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

        #t_test = (L_F*M_TBT*rho_p*V_p)/(r_TBT * M_unit)

        dlpdt = (M_unit * r_TBT) / (M_TBT * rho_p * V_p)

        dlcdt = (M_Cu2O*r_Cu2O)/(2*V_c*rho_c) 

        return np.array([dlpdt, dlcdt])

        
# Run method
fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)
tol = 10e-20
n = 100
Nmax = 10**(n+10)


### 400 days, Bueno Aires coastline
def p_400_temp(t):
    return temperatures[0]

def p_400_sal(t):
    return salinity[0]

def p_400_pH(t):
    return pH[0]

model_400 = ADVANCED_MODEL(p_400_temp, p_400_sal, p_400_pH)

t0_400, T_end_400 = 0, 400

l0_400 = np.array([0, 10**(-n)])

t_400, y_400 = fehlberg(l0_400, t0_400, T_end_400, model_400, Nmax, tol)


def plot_conversion_400():
    plt.plot(t_400, y_400/L_F)
    plt.legend(["$X_p$", "$X_c$"])
    plt.grid(True)
    plt.ylabel("Conversion [-]")
    plt.xlabel("Time [days]")
    plt.title("Conversion vs time, advanced model. Normal model used.")
    plt.show()


def plot_length_400():
    plt.plot(t_400, y_400*10**3)
    plt.legend(["$l_p$", "$l_c$"])
    plt.grid(True)
    plt.ylabel("Length of the moving fronts [mm]")
    plt.xlabel("Time [days]")
    plt.title("Length of front, advanced model. Normal model used.")
    plt.show()


def thickness_leached_layer_400():
    thickness = np.zeros(len(t_400))
    for i in range(len(t_400)):
        thickness[i] = y_400[i][1] - y_400[i][0]

    plt.plot(t_400, thickness*10**3)
    plt.xlabel("Time [days]")
    plt.ylabel("Thickness [mm]")
    plt.title("Thickness of the leached layer, advanced. Normal model used.")
    plt.grid(True)
    plt.show()

# plot_conversion_400()
# # plot_length_400()
# thickness_leached_layer_400()


### THE NEXT 20 DAYS
model_20 = ADVANCED_MODEL(p_temperature, p_salinity, p_pH)

t0_20, T_end_20 = 0, 19.9

l0_20 = np.array([y_400[-1][0], y_400[-1][1]])

t_20, y_20 = fehlberg(l0_20, t0_20, T_end_20, model_20, Nmax, tol)


def plot_conversion_20():
    plt.plot(t_20, y_20/L_F)
    plt.legend(["$X_p$", "$X_c$"])
    plt.grid(True)
    plt.ylabel("Conversion [-]")
    plt.xlabel("Time [days]")
    plt.title("Conversion vs time, advanced model. Normal model used.")
    plt.show()

def plot_length_20():
    plt.plot(t_20, y_20*10**3)
    plt.legend(["$l_p$", "$l_c$"])
    plt.grid(True)
    plt.ylabel("Length of moving fronts [mm]")
    plt.xlabel("Time [days]")
    plt.title("Length of moving fronts, voyage, advanced model. Normal model used.")
    plt.show()

def plot_thickness_20():
    thickness = np.zeros(len(t_20))
    for i in range(len(t_20)):
        thickness[i] = y_20[i][1] - y_20[i][0]

    plt.plot(t_20, thickness*10**3)
    plt.xlabel("Time [days]")
    plt.ylabel("Thickness [mm]")
    plt.title("Thickness of the leached layer, voyage, advanced. Normal model used.")
    plt.grid(True)
    plt.show()


# plot_conversion_20()
# # plot_length_20()
# plot_thickness_20()

### 400 + 20 days modeling
t_1 = t_400[:]
y_1 = y_400[:]

t_2 = t_20[:] + 400
y_2 = y_20[:]

t_3 = np.concatenate((t_1, t_2))
y_3 = np.concatenate((y_1, y_2))

plt.plot(t_3, y_3/L_F)
plt.legend(["$X_p$", "$X_c$"])
plt.xlabel("Time [days]")
plt.ylabel("Conversion [-]")
plt.grid(True)
plt.show()



