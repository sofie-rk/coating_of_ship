import numpy as np
from numpy.linalg import norm

class EmbeddedExplicitRungeKutta:
    def __init__(self, a, b, c,  bhat=None, order=None):
        # Butcher Table
        self.a = a
        self.b = b
        self.c = c
        self.bhat = bhat
        self.order = order


    def __call__(self, y0, t0, T, f, Nmax, tol=1e-3):
        # Extract Butcher table
        a, b, c, bhat, order = self.a, self.b, self.c, self.bhat, self.order

        # Some parameters controlling the time-step choice (Machine precision)
        fac = 0.8
        facmax = 2
        facmin = 0.1
        err  = 0
        
        # Stages (RK method)
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

            dyhat = np.zeros_like(y, dtype=np.double)
            for j in range(s):
                dyhat += bhat[j]*ks[j]

            # Error estimate
            err = dt*norm(dy - dyhat)

            # Accept time-step
            if err <= tol:
                y_next = y + dt*dyhat
                t_next = t + dt

                if (y_next[0] > 1):
                    y_next[0] = 1
                
                if (y_next[1] > 1):
                    y_next[1] = 1

                # if (t_next > 1):
                #     t_next = 1

                ys.append(y_next)
                ts.append(t_next)    
                
            else: # rejected
                print(f"Step is rejected at t = {t} with err = {err}")
                N_rej += 1
                ys_rej.append(y + dt*dyhat)
                ts_rej.append(t + dt)
            

            dt = min(dt*min(facmax, max(facmin, fac*(tol/err)**(1/(order)))),abs(T-t))
            

        
        print(f"Finishing time-stepping reaching t = {ts[-1]} with final time T = {T}")
        print(f"Used {N} steps out of {Nmax} with {N_rej} being rejected")
          
        
        return (np.array(ts), np.array(ys))

# end of class EmbeddedExplicitRungeKutta

### FEHLBERG-RUNGEKUTTA BUTCHER TABLEAU
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