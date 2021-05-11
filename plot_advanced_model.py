from advanced_model import t_400, y_400, t_20, y_20
from labels import *

### FIRST 400 DAYS ###

def plot_conversion_400():
    plt.plot(t_400, y_400/L_F)
    plt.legend(["$X_p$", "$X_c$"])
    plt.grid(True)
    plt.ylabel("Conversion [-]")
    plt.xlabel("Time [days]")
    plt.title("Conversion vs time, advanced model. Normal model used.")
    plt.show()


def plot_length_400():
    plt.plot(t_400, y_400*10**3)
    plt.legend(["$l_p$", "$l_c$"])
    plt.grid(True)
    plt.ylabel("Length of the moving fronts [mm]")
    plt.xlabel("Time [days]")
    plt.title("Length of front, advanced model. Normal model used.")
    plt.show()


def thickness_leached_layer_400():
    thickness = np.zeros(len(t_400))
    for i in range(len(t_400)):
        thickness[i] = y_400[i][1] - y_400[i][0]

    plt.plot(t_400, thickness*10**3)
    plt.xlabel("Time [days]")
    plt.ylabel("Thickness [mm]")
    plt.title("Thickness of the leached layer, advanced. Normal model used.")
    plt.grid(True)
    plt.show()

# plot_conversion_400()
# # plot_length_400()
# thickness_leached_layer_400()


### 20 day voyage

