from values import *
import scipy.integrate
#from scipy.integrate import ode

def model(t, z):

    zp = z[0]
    zc = z[1]

    dzp_dt = 1

    dzc_dt = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/(zc*L_F - zp*L_F)

    return np.array([dzp_dt, dzc_dt])

a = scipy.integrate.ode(model).set_integrator('vode', method='bdf', order=15)


t0 = 0
y0 = [0, 10**(-10)]
a.set_initial_value(y0,t0)

a.successful(self)

print(a)