import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from values import*

class EmbeddedExplicitRungeKutta:
    def __init__(self, a, b, c,  bhat=None, order=None):
        self.a = a
        self.b = b
        self.c = c
        self.bhat = bhat
        self.order = order


    def __call__(self, y0, t0, T, f, Nmax, tol=1e-3):
        # Extract Butcher table
        a, b, c, bhat, order = self.a, self.b, self.c, self.bhat, self.order

        # Some parameters controlling the time-step choice 
        # Machine precision
        eps = 1e-15
        fac = 0.8
        facmax = 5.0
        facmin = 0.1
        err  = 0
        
        # Stages
        s = len(b)
        ks = [np.zeros_like(y0, dtype=np.double) for s in range(s)]

        # Start time-stepping
        ys = [y0]
        ts = [t0]
        
        # Store rejected time-steps
        ts_rej = []
        ys_rej = []
        dt = (T - t0)/Nmax
        # Counting steps
        N = 0
        N_rej = 0
        
        while(ts[-1] < T and N < Nmax):
            t, y = ts[-1], ys[-1]
            N += 1
            
#             print("---------------------------------------------------")
#             print(f"Compute tentative step N = {N} starting from t = {t} with dt = {dt}")
            
            # Compute stages derivatives k_j
            for j in range(s):
                t_j = t + c[j]*dt
                dY_j = np.zeros_like(y, dtype=np.double)
                for l in range(j):
                    dY_j += a[j,l]*ks[l]

                ks[j] = f(t_j, y + dt*dY_j)
                
            # Compute next time-step
            dy = np.zeros_like(y, dtype=np.double)
            for j in range(s):
                dy += b[j]*ks[j]

            if bhat is None:
                ys.append(y + dt*dy)
                ts.append(t + dt)
            else:
                dyhat = np.zeros_like(y, dtype=np.double)
                for j in range(s):
                    dyhat += bhat[j]*ks[j]

                # Error estimate
#                 err = max(dt*norm(dy - dyhat), norm(y)*eps)
                err = dt*norm(dy - dyhat)

                # Accept time-step
#                 if True:
                if err <= tol:
                    ys.append(y + dt*dyhat)
                    ts.append(t + dt)
                else:
                    print(f"Step is rejected at t = {t} with err = {err}")
                    N_rej += 1
                    ys_rej.append(y + dt*dyhat)
                    ts_rej.append(t + dt)
                
                # New step size
                dt = min(dt*min(facmax, max(facmin, fac*(tol/err)**(1/(order)))),abs(T-t))
#                 dt = 0.8*(tol/err)**(1/(order))*dt
#                 print(f"New dt = {dt}")
#                 print(f"Error is err = {err}")
        
        print(f"Finishing time-stepping reaching t = {ts[-1]} with final time T = {T}")
        print(f"Used {N} steps out of {Nmax} with {N_rej} being rejected")
          
        
        return (np.array(ts), np.array(ys))

# end of class EmbeddedExplicitRungeKutta


### FEHLBERG-RUNGEKUTTA DEF
a =  np.array([[0.0, 0,   0,   0,   0],
                       [1/2, 0,   0,   0,   0],
                       [0, 1/2,   0,   0,   0],
                       [0,   0,   1,   0,   0],
                       [1/6, 1/3, 1/3, 1/6, 0],
                        ])
                    
b = np.array([1/6, 1/3, 1/3, 0, 1/6])
bhat = np.array([1/6, 1/3, 1/3, 1/6, 0])

c = np.array([0.,
                       1/2,
                       1/2,
                       1,
                       1])

order =  4



class MODEL:

    # def __init__(self, M_TBT):
    #     self.M_TBT = M_TBT

    def __call__(self, t, z):
        dzp_dt = 1
        dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(z[1]*L_F - z[0]*L_F)

        return np.array([dzp_dt, dzc_dt])


# Define a model
model = MODEL()

n = 100
z0 = np.array([0, 10**(-n)])


# Run method

# Initial time, stop time
t0, T = 0, 1
Nmax = 10**(n-1)

# Tolerance in adaptive method
tol = 1.0e-3


fehlberg = EmbeddedExplicitRungeKutta(a, b, c, bhat, order)


ts, ys = fehlberg(z0, t0, T, model, Nmax, tol)

plt.plot(ts, ys, '-o')
plt.legend(["$z_p$", "$z_c$"])
plt.xlabel("Tau")
plt.show()

