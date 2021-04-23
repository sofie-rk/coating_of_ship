from values import *

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


print("STARTING ODE SOLVER")

def model(z, tau):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (rho_p*V_p*M_TBT)/(r_TBT*M_unit) * D_CuCl/(zc*L_F - zp*L_F) * 1/(2*V_c) * alpha

    return [dzp_dt, dzc_dt]

def ODEsolver():

    # Initial conditions
    init = [0, 0]

    n = 100

    tau_list = np.linspace(0, 1, n)

    # Store solution
    zp_list = np.zeros(len(tau_list))
    zc_list = np.zeros(len(tau_list))

    # Add inital conditions
    zp_list[0] = init[0]
    zc_list[1] = init[1]

    for i in range(1, n):
        tspan = [tau_list[i-1], tau_list[i]]

        z  = odeint(model, init, tspan)

        zp_list[i] = z[1][0]
        zc_list[i] = z[1][1]

        # New init conditions
        init = z[1]

    plt.plot(tau_list, zp_list, label="$z_p$")
    plt.plot(tau_list, zc_list, label="$z_c$")
    plt.xlabel("Time t")
    plt.legend()
    plt.show()

    #print(zp_list)
    print(zc_list)

ODEsolver()





