from conversion_vs_tau import ys, ts
from values import L_F

import matplotlib.pyplot as plt
import numpy as np



thickness = np.zeros(len(ts))
for i in range(len(ts)):
    thickness[i] = (ys[i][1] - ys[i][0]) * L_F * 10**3  #[mm]
plt.plot(ts, thickness, '.')
plt.title("Thickness of leached layer ")
plt.xlabel("Tau")
plt.ylabel("$l_c - l_p$ [mm]")
plt.show()
