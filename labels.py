x_label_day = "Time [days]"
tau_label = "Tau [-]"
conversion_label = "Conversion [-]"
conversion_legend = ["$X_p$", "$X_c$"]
release_Cu_label = r"Release of $Cu^{2+}$  [$\mu$ g / $cm^2$ day]"
Xp_label = "$X_p$"
Xc_label = "$X_c$"
thickness_label = "Thickness of leached layer [mm]"

def temp_label(T):
    return "T = " + str(T-273) + "$^o$C"

def salinity_label(sal):
    return "Salinity = " + str(sal) + "g salt / kg seawater"

def pH_label(pH):
    return "pH = " + str(pH)