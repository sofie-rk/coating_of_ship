from values import *

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


print("STARTING ODE SOLVER")

def model(z, tau):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (M_TBT*M_Cu2O/M_unit) * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * (C_CuCl_s-0)/(zc*L_F-zp*L_F)

    return [dzp_dt, dzc_dt]

def ODEsolver():

    # Initial conditions
    init = [0, 10**(-15)]

    n = 51

    tau_list = np.linspace(0, 1, n)

    # Store solution
    zp_list = np.zeros(len(tau_list))
    zc_list = np.zeros(len(tau_list))

    # Add inital conditions
    zp_list[0] = init[0]
    zc_list[0] = init[1]

    for i in range(1, n):
        tspan = [tau_list[i-1], tau_list[i]]

        z  = odeint(model, init, tspan)

        zp_list[i] = z[1][0]
        zc_list[i] = z[1][1]

        # New init conditions
        init = z[1]

    plt.plot(tau_list, zp_list, label="$z_p$")
    plt.plot(tau_list, zc_list, label="$z_c$")
    plt.xlabel("Tau")
    plt.legend()
    plt.show()


ODEsolver()





