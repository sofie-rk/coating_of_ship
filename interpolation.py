import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d


time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

temperatures = [15.94+273, 17.72+273, 19.38+273, 20.92+273, 22.32+273, 23.56+273, 24.64+273, 25.53+273, 26.21+273, 26.68+273, 26.92+273, 26.91+273, 26.62+273, 26.07+273, 25.21+273, 24.04+273, 22.55+273, 20.71+273, 18.51+273, 15.93+273, 12.96+273]

salinity = [35.46, 35.38, 35.57, 35.76, 35.81, 35.69, 35.45, 35.15, 34.86, 34.67, 34.61, 34.70, 34.92, 35.21, 35.49, 35.67, 35.67, 35.43, 34.96, 34.37, 33.87]

pH = [8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 7.95, 7.78, 7.64, 7.53, 7.4] 

p_temperature = interp1d(time, temperatures)
p_salinity = interp1d(time, salinity)
p_pH = interp1d(time, pH)

# plt.plot(time, p_temperature(time), label="interpolation")
# plt.plot(time, temperatures, 'o')
# plt.legend()
# plt.xlabel("Time [days]")
# plt.ylabel("Temperature [K]")
# plt.show()
