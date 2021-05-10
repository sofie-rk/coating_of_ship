from simple_model import tau_simple, X_simple
from simple_model import tf_simple
from values import L_F
from values import D_CuCl, C_CuCl_s, M_Cu

from labels import x_label_day

import numpy as np
import matplotlib.pyplot as plt

def plot_conversion_simple():
    plt.plot(tau_simple*tf_simple, X_simple)
    plt.legend(["$X_p$", "$X_c$"])
    plt.ylabel("Conversion [-]")
    plt.xlabel(x_label_day)
    plt.title("Simple model, conversion vs time. Dimensionless model used.")
    plt.grid(True)
    plt.show()


### THICKNESS OF THE LEACHED LAYER ###

def plot_thickness_simple():
    thickness = np.zeros(len(tau_simple))
    for i in range(len(tau_simple)):
        thickness[i] = (X_simple[i][1] - X_simple[i][0]) * L_F * 10**3  #[mm]
    plt.plot(tau_simple[:]*tf_simple, thickness[:])
    plt.title("Thickness of leached layer. Simple model. Dimensionless model used.")
    plt.xlabel(x_label_day)
    plt.ylabel("$l_c - l_p$ [mm]")
    plt.show()

plot_conversion_simple()
plot_thickness_simple()


### RELEASE OF Cu2+

def plot_release_Cu():

    r_Cu2O_simple = np.zeros(len(tau_simple))

    for i in range(len(r_Cu2O_simple)):

        r_Cu2O_simple[i] = 2 * D_CuCl*C_CuCl_s/(L_F*(X_simple[i][1] - X_simple[i][0])) * M_Cu * 10**(-4) # [micro g/cm2 day]
    
    plt.plot(tau_simple[1500:-10]*tf_simple, r_Cu2O_simple[1500:-10])
    plt.xlabel(x_label_day)
    plt.ylabel(r"Release of $Cu^{2+}$  [$\mu$ g / $cm^2$ day]")
    plt.ylim(20, 40)
    plt.show()



plot_release_Cu()