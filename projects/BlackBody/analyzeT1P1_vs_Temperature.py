"""
Used to analyze the up/down transition rate based on the T1, P1 as
a function of blackbody temperature
"""

import numpy as np
import matplotlib.pyplot as plt

# array for each temp in mK, P1 in units of %, T1 in units of us, PSD in units of ms
units = [1e-3, 0.01, 1e-6]
T1_other = 24.0 * 1e-6
data_2d = np.asarray([
    [76.0, 1.629, 20.75], [82.5, 1.790, 21.30], [95.5, 1.768, 20.50],
    [104.8, 1.838, 20.66], [110.6, 2.098, 21.93], [116.5, 1.976, 23.86],
    [125, 2.138, 21.21], [138.5, 2.493, 20.25], [145.5, 2.614],
    [155, 2.602], [163, 2.899], [170, 3.294],
    [175, 3.550], [180, 3.736], [186, ]
])

for i, d in enumerate(data_2d):  # convert units to SI
    data_2d[i] = np.multiply(d, units)

# make data array
Temp_array = data_2d[:, 0]
P1_array = data_2d[:, 1]
T1_total_array = data_2d[:, 2]

# extract QP induced part
T1_QP_array = [1 / (1 / T - 1 / T1_other) for T in T1_total_array]
Up_total_array = T1_total_array * (1 / P1_array - 1)  # Upward transition rate

# convert time to rates
Gamma_down_total = [1 / T for T in T1_total_array]
Gamma_down_QP = [1 / T for T in T1_QP_array]
Gamma_up_QP = [1 / T for T in Up_total_array]

Gamma_down_QP_scale = [np.log(1 / G) for G in Gamma_down_QP]
Gamma_up_QP_scale = [np.log(1 / G) for G in Gamma_up_QP]

# plt.plot(Temp_array, P1_array, label='P1')
# plt.plot(Temp_array, Up_total_array, label='Up')
# plt.plot(Temp_array, T1_QP_array, label='T1_QP')
# plt.plot(Temp_array, T1_total_array, label='T1_tot')
# plt.yscale("log")
# plt.xlabel('Temp (Kelvin)')
# plt.ylabel('Time (sec)')
# plt.grid()
# plt.legend()
# plt.show()

"""Plot rate"""
plt.plot(Temp_array, Gamma_up_QP, label='1/Up')
plt.plot(Temp_array, Gamma_down_QP, label='1/T1_QP')
plt.plot(Temp_array, Gamma_down_total, label='1/T1_tot')
plt.yscale("log")
plt.xlabel('Temp (Kelvin)')
plt.ylabel('Rate (1/sec)')
plt.grid()
plt.legend()
plt.show()


"""Processed data and function fit"""
# plt.scatter([1/T for T in Temp_array], Gamma_down_QP_scale, label='QP_down')
# plt.scatter([1 / T for T in Temp_array], Gamma_up_QP_scale, label='QP_Up')
# plt.scatter([1 / T for T in Temp_array], [0.5 / T - 10.5 for T in Temp_array],
#             label='QP_fit')
# # plt.scatter([1/T for T in Temp_array[:27]], Gamma_PSD_scale[:27], label='PSD')
# plt.xlabel('1/Temp(Kelvin)')
# plt.ylabel('ln(df/Gamma)')
# plt.grid()
# plt.legend()
# plt.show()
