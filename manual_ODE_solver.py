from values import *

import numpy as np
import matplotlib.pyplot as plt

def model(t, z):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(zc*L_F - zp*L_F)

    return np.array([dzp_dt, dzc_dt])

def explicit_euler(y0, t0, T, f, Nmax):
    y_sol = [y0] # initial condition
    t_sol = [t0]

    h = (T-t0)/Nmax # step length

    while (t_sol[-1] < T):

        t, y = t_sol[-1], y_sol[-1] 

        y_sol.append(y + h*f(t, y))
        t_sol.append(t + h)

    return np.array(t_sol), np.array(y_sol)

def heun(y0, t0, T, f, Nmax):
    y_sol = [y0]
    t_sol = [t0]
    h = (T-t0)/Nmax
    
    while(t_sol[-1] < T):
        t, y = t_sol[-1], y_sol[-1]

        k1 = f(t, y)
        k2 = f(t+h, y+h*k1)

        y_sol.append(y + 0.5*h*(k1+k2))
        t_sol.append(t+h)

    return np.array(t_sol), np.array(y_sol)

def RK4(y0, t0, T, f, Nmax):

    # similar to ode45 in matlab

    y_sol = [y0]
    t_sol = [t0]
    h = (T-t0)/Nmax
    
    while(t_sol[-1] < T):
        t,y = t_sol[-1], y_sol[-1]

        k1 = f(t, y)
        k2 = f(t + h/2, y + h*k1/2)
        k3 = f(t + h/2, y + h*k2/2)
        k4 = f(t + h/2, y + h*k3)

        y_sol.append(y + 1/6*h *(k1 + 2*k2 + 2*k3 + k4))
        t_sol.append(t + h)

    return np.array(t_sol), np.array(y_sol)


tau_0, tau_T = 0, 1
z0 = np.array([0, 10**(-5)]) # zp_0 and zc_0

Nmax = 5000

tau_sol, z_sol = RK4(z0, tau_0, tau_T, model, Nmax)


plt.plot(tau_sol, z_sol, label="$z_p$")
plt.title("zc(tau=0) = " + str(z0[-1]))
plt.legend(["$z_p$", "$z_c$"])
plt.xlabel("Tau")
plt.show()