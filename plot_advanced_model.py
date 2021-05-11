from advanced_model import t_400, y_400, t_20, y_20, t_420, y_420
from advanced_model import r_Cu2O_advanced
from conditions_advanced import *
from labels import *
from values import L_F, M_Cu

import numpy as np
import matplotlib.pyplot as plt

### FIRST 400 DAYS ###

def plot_conversion_400():
    plt.plot(t_400, y_400/L_F)
    plt.legend(["$X_p$", "$X_c$"])
    plt.grid(True)
    plt.ylabel("Conversion [-]")
    plt.xlabel("Time [days]")
    plt.show()

    #print("Advanced, conversion after 400 days: ", y_400[-1]/L_F)


def plot_length_400():
    plt.plot(t_400, y_400*10**3)
    plt.legend(["$l_p$", "$l_c$"])
    plt.grid(True)
    plt.ylabel("Length of the moving fronts [mm]")
    plt.xlabel("Time [days]")
    plt.show()


def thickness_leached_layer_400():
    thickness = np.zeros(len(t_400))
    for i in range(len(t_400)):
        thickness[i] = y_400[i][1] - y_400[i][0]

    plt.plot(t_400, thickness*10**3)
    plt.xlabel(x_label_day)
    plt.ylabel("Thickness of leached layer [mm]")
    plt.grid(True)
    plt.show()

def release_Cu_400():
    r_Cu = np.zeros(len(t_400))

    for i in range(len(r_Cu)):

        t = t_400[i]

        r_Cu[i] = r_Cu2O_advanced(p_400_temp(t), p_400_pH(t), p_400_sal(t), y_400[i][1], y_400[i][0]) * M_Cu * 10**(-4) # [micro g/cm2 day]

    
    plt.plot(t_400, r_Cu)
    plt.xlabel(x_label_day)
    plt.ylabel(release_Cu_label)
    plt.grid(True)
    plt.ylim(0,60)
    plt.show()

    #print("Advanced, release rate after 400 days: ", r_Cu[-1])


### 20 day voyage

def plot_conversion_20():
    plt.plot(t_20, y_20/L_F)
    plt.legend(conversion_legend)
    plt.grid(True)
    plt.ylabel(conversion_label)
    plt.xlabel(x_label_day)
    plt.show()

def plot_length_20():
    plt.plot(t_20, y_20*10**3)
    plt.legend(["$l_p$", "$l_c$"])
    plt.grid(True)
    plt.ylabel("Length of moving fronts [mm]")
    plt.xlabel(x_label_day)
    plt.title("Length of moving fronts, voyage, advanced model. Normal model used.")
    plt.show()

def plot_thickness_20():
    thickness = np.zeros(len(t_20))
    for i in range(len(t_20)):
        thickness[i] = y_20[i][1] - y_20[i][0]

    plt.plot(t_20, thickness*10**3)
    plt.xlabel(x_label_day)
    plt.ylabel("Thickness [mm]")
    plt.grid(True)
    plt.show()

def release_Cu_20():
    r_Cu = np.zeros(len(t_20))

    for i in range(len(r_Cu)):

        t = t_20[i]

        r_Cu[i] = r_Cu2O_advanced(p_20_temp(t), p_20_pH(t), p_20_sal(t), y_20[i][1], y_20[i][0]) * M_Cu * 10**(-4) # [micro g/cm2 day]

    
    plt.plot(t_20, r_Cu)
    plt.xlabel(x_label_day)
    plt.ylabel(release_Cu_label)
    plt.grid(True)
    plt.show()


### 420 days

def plot_conversion_420():

    plt.plot(t_420, y_420/L_F)
    plt.legend(conversion_legend)
    plt.xlabel(x_label_day)
    plt.ylabel(conversion_label)
    plt.grid(True)
    plt.show()

def release_Cu_420():
    r_Cu = np.zeros(len(t_420))

    for i in range(len(t_420)):

        t = t_420[i]

        if t<400:
            r_Cu[i] = r_Cu2O_advanced(p_400_temp(t), p_400_pH(t), p_400_sal(t), y_420[i][1], y_420[i][0]) * M_Cu * 10**(-4) # [micro g/cm2 day]
        else:
            r_Cu[i] = r_Cu2O_advanced(p_20_temp(t-400), p_20_pH(t-400), p_20_sal(t-400), y_420[i][1], y_420[i][0]) * M_Cu * 10**(-4) # [micro g/cm2 day]

    
    plt.plot(t_420, r_Cu)
    plt.xlabel(x_label_day)
    plt.ylabel(release_Cu_label)
    plt.grid(True)
    plt.ylim(0,60)
    plt.show()

    print("Release after 420 days: ", r_Cu[-1])




plot_conversion_400()
# # plot_length_400()
thickness_leached_layer_400()
release_Cu_400()

plot_conversion_20()
# # plot_length_20()
plot_thickness_20()
release_Cu_20()

plot_conversion_420()
release_Cu_420()