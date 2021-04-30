def explicit_euler(y0, t0, T, f, Nmax):
    y_sol = [y0] # initial condition
    t_sol = [t0]

    h = (T-t0)/Nmax # step length

    while (t_sol[-1] < T):

        t, y = t_sol[-1], y_sol[-1] 

        y_sol.append(y + h*f(t, y))
        t_sol.append(t + h)

    return np.array(t_sol), np.array(y_sol)

def heun(y0, t0, T, f, Nmax, tol):
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
    # Standard RK4 code

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
