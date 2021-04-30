rho_p = 1200            # [kg/m3] 
rho_c = 6000            # [kg/m3]
L_F = 0.55 * 10**(-3)   # [m]

M_Cu2O = 0.143          # [kg Cu2O/mol]
M_unit = 0.574          # [kg unit/mol] 
M_TBT = 0.2897          # [kg TBT/mol]  

r_TBT = 4.2 * 10**(-9) * 10**4  # [kg TBT/m^2 day]

V_p = 0.7   # volume fraction of polymer
V_c = 0.3   # volume fraction of pigment

D_CuCl = 2.0 * 10**(-12) * 3600 * 24    # [m2/day] diffusivity of CuCl2 through the leached layer

C_CuCl_s = 8.6 * 10**(-2)   # [mol/m3] concentration of CuCl- at polymer/leached layer interface


tf = L_F*M_TBT*rho_p*V_p/(r_TBT*M_unit)


const = (M_TBT*M_Cu2O)/M_unit * (rho_p*V_p*D_CuCl)/(r_TBT*2*V_c*rho_c) * C_CuCl_s/L_F

