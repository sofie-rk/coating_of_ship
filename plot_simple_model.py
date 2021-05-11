from simple_model import tau_simple, X_simple
from simple_model import tf_simple
from values import L_F
from values import D_CuCl, C_CuCl_s, M_Cu

from labels import *

import numpy as np
import matplotlib.pyplot as plt

def plot_conversion_simple():
    # Plotting conversion of polymer can Cu2O vs time
    plt.plot(tau_simple*tf_simple, X_simple)
    plt.legend(conversion_legend)
    plt.ylabel(conversion_label)
    plt.xlabel(x_label_day)
    plt.grid(True)
    plt.show()


def plot_thickness_simple():
    # Plotting thickness of leached layer vs time
    thickness = np.zeros(len(tau_simple))
    for i in range(len(tau_simple)):
        thickness[i] = (X_simple[i][1] - X_simple[i][0]) * L_F * 10**3  #[mm]
    plt.plot(tau_simple[:]*tf_simple, thickness[:])
    plt.title("Thickness of leached layer. Simple model. Dimensionless model used.")
    plt.xlabel(x_label_day)
    plt.ylabel("$l_c - l_p$ [mm]")
    plt.show()


### RELEASE OF Cu2+

def plot_release_Cu():
    r_Cu2O_simple = np.zeros(len(tau_simple))

    for i in range(len(r_Cu2O_simple)):

        r_Cu2O_simple[i] = D_CuCl*C_CuCl_s/(L_F*(X_simple[i][1] - X_simple[i][0])) * M_Cu * 10**(-4) # [micro g/cm2 day]
    
    plt.plot(tau_simple[:]*tf_simple, r_Cu2O_simple[:])
    plt.xlabel(x_label_day)
    plt.ylabel(release_Cu_label)
    plt.ylim(0, 50)
    plt.grid(True)
    plt.show()

plot_conversion_simple()
plot_release_Cu()
