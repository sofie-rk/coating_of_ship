from values import *

import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import norm

def model(t, z):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(zc*L_F - zp*L_F)

    return np.array([dzp_dt, dzc_dt])

def heun_euler(y0, t0, T, f, Nmax, tol=1e-3):

    # Store inital condition
    ys = [y0]
    ts = [t0]

    # step size 
    h = (T-t0)/Nmax

    # number of steps
    s = 1 # (euler)

    order = 1 # ???
    
    # Counting
    N = 0
    N_rejected = 0

    while (ts[-1] < T and N < Nmax):
        t, y = ts[-1], ys[-1]   # access previous solition 
        N += 1  

        # Parameters controlling time-step choice
        eps = 1e-15 # ???
        fac = 0.8
        facmax = 5
        facmin = 0.1
        err = 0


        # Calculate k1 and k2
        k1 = f(t, y)
        k2 = f(t+h, y+ h*k1)

        # Compute next step
        y_next_1 = y + h*k1
        y_next_2 = y + h/2*(k1 + k2)

        # Error estimate
        err = h*norm(y_next_1 - y_next_2)

        # Accept time step
        if (err <= tol):
            ys.append(y + h*y_next_2)
            ys.append(t + h)

        else:
            print(f"Step is rejected at t = {t} with error = {err}")

            N_rej += 1

            # new step size
            h = 0.8*(tol/err)**(1/order)*h
            print("New h: ", h)

        return np.array(ts), np.array(ys)



tau_0, tau_T = 0, 1
z0 = np.array([0, 10**(-7)]) # zp_0 and zc_0

Nmax = 100000

tau_sol, z_sol = heun_euler(z0, tau_0, tau_T, model, Nmax)




