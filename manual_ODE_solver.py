from values import *

import numpy as np
import matplotlib.pyplot as plt
import timeit

def model(t, z):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(zc*L_F - zp*L_F)

    return np.array([dzp_dt, dzc_dt])


def RK4(y0, t0, T, f, Nmax):
    # Standard RK4 code

    t_sol = [t0]
    y_sol = [y0]
    
    h = (T-t0)/Nmax
    
    while(t_sol[-1] < T):
        t,y = t_sol[-1], y_sol[-1]

        k1 = f(t, y)
        k2 = f(t + h/2, y + h*k1/2)
        k3 = f(t + h/2, y + h*k2/2)
        k4 = f(t + h, y + h*k3)

        y_sol.append(y + 1/6*h *(k1 + 2*k2 + 2*k3 + k4))
        t_sol.append(t + h)

    return np.array(t_sol), np.array(y_sol)

n = 8

tau_0, tau_T = 0, 1
z0 = np.array([0, 10**(-n)]) # zp_0 and zc_0

Nmax = 10**(n-1)


print("STARTING PROGRAM")
start_time = timeit.default_timer()
tau_sol, z_sol = RK4(z0, tau_0, tau_T, model, Nmax)
end_time = timeit.default_timer()
print(f"FINISHED USING t = {end_time-start_time} seconds")


plt.plot(tau_sol, z_sol, label="$z_p$")
plt.title(f"zc(tau=0) = {z0[-1]}  Nmax = {Nmax}")
plt.legend(["$z_p$", "$z_c$"])
plt.xlabel("Tau")
plt.show()