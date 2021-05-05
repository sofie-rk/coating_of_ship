from numpy import exp

seconds_in_a_day = 3600 * 24

rho_p = 1200            # [kg/m3] 
rho_c = 6000            # [kg/m3]
L_F = 0.55 * 10**(-3)   # [m]

M_Cu2O = 0.143          # [kg Cu2O/mol]
M_unit = 0.574          # [kg unit/mol] 
M_TBT = 0.2897          # [kg TBT/mol]  

r_TBT = 4.2 * 10**(-9) * 10**4  # [kg TBT/m^2 day]

V_p = 0.7   # volume fraction of polymer
V_c = 0.3   # volume fraction of pigment

D_CuCl = 2.0 * 10**(-12) * seconds_in_a_day    # [m2/day] diffusivity of CuCl2 through the leached layer

C_CuCl_s = 8.6 * 10**(-2)   # [mol/m3] concentration of CuCl- at polymer/leached layer interface


tf = L_F*M_TBT*rho_p*V_p/(r_TBT*M_unit)


const = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/L_F



### FOR PART 6

L_ship = 100    # [m]
mu_seawater = 9 * 10**(-4)     # [kg/ms]
rho_seawater = 1025             # [kg/m3]
u_ship = 15 * 0.514444          # [m/s]

D_CuCl_sw = 0.9 * 10**(-9) * seconds_in_a_day     # [m2/day]

Re = L_ship * u_ship * rho_seawater / mu_seawater # Reynolds number
print(Re)

k_L = 0.0365*Re**(0.8) * D_CuCl_sw / L_ship

print("k_L: ", k_L)

def k1(T):
    k0 = 1.63 * seconds_in_a_day           # [mol (mol/m3)**(-1.32) / m2 day]
    Ea = 61.1*1000   # [J/mol]
    R = 8.314        # [J/Kmol]

    return k0 * exp(-Ea / (R*T))

k2 = 8      # [(mol/m3)**(-0.43)]

def k_Cu2O(T):
    k0 = 146.6 * seconds_in_a_day     # [mol (mol/m3)**(-3) / m2 day]
    Ea = 50.2 * 1000    # [J/mol]
    R = 8.314           # [J/Kmol]

    return k0 * exp(-Ea/(R*T))

k_1 = 5.1*10**(-6) * seconds_in_a_day      # [m/day]

