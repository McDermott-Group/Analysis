"""
Used to analyze the up/down transition rate based on the T1, P1 as
a function of blackbody temperature
"""

from antennalib import getGamma_pa

import numpy as np
import matplotlib.pyplot as plt

# array for each temp in mK, P1 in units of %, T1 in units of us, PSD in units of ms
units = [1e-3, 0.01, 1e-6, 1e-3]
T1_other = 25.0*1e-6
# Temp_P1_T1_1 = [80, 2.1, 20.8, 3.145]
# Temp_P1_T1_2 = [185, 4.3, 17.6, 2.847]
# Temp_P1_T1_3 = [200, 4.9, 16.1, 2.750]
# Temp_P1_T1_4 = [210, 5.6, 19.2, 2.709]
# Temp_P1_T1_5 = [220, 6.0, 19.0, 2.753]
# Temp_P1_T1_6 = [229, 6.2, 15.5, 2.683]
# Temp_P1_T1_7 = [237, 6.9, 16.8, 2.608]
# Temp_P1_T1_8 = [245, 7.1, 18.9, 2.649]
# Temp_P1_T1_9 = [252, 7.8, 18.2, 2.602]
# Temp_P1_T1_10 = [260, 8.1, 16.0, 2.441] # below is taken at different days
# Temp_P1_T1_11 = [266, 8.2, 16.5, 2.437]
# Temp_P1_T1_12 = [274, 8.4, 16.2, 2.469]
# Temp_P1_T1_13 = [284, 8.9, 17.0, 2.383]
# Temp_P1_T1_14 = [294, 9.3, 14.2, 2.228]
# Temp_P1_T1_15 = [305, 9.9, 16.6, 2.228] #question
# Temp_P1_T1_16 = [318, 10.4, 17.6, 2.165]
# Temp_P1_T1_17 = [331, 11.3, 15.1, 1.993]
# Temp_P1_T1_18 = [342, 11.4, 14.4, 1.776]
# Temp_P1_T1_19 = [353, 12.2, 12.2, 1.526]
# Temp_P1_T1_20 = [364, 12.2, 14.9, 1.249]
# Temp_P1_T1_21 = [375, 13.7, 13.3, 1.049]
# Temp_P1_T1_22 = [386, 13.3, 12.2, 0.810]    # PSD above not so much trust   100us sampling rate
# Temp_P1_T1_23 = [396, 14.9, 13.5, 0.654]
# Temp_P1_T1_24 = [404, 15.3, 13.4, 0.506]
# Temp_P1_T1_25 = [412, 15.4, 12.8, 0.422]
# Temp_P1_T1_26 = [424, 17.0, 11.7, 0.355]
# Temp_P1_T1_27 = [438, 17.2, 11.4, 0.3]
# Temp_P1_T1_28 = [449, 19.4, 10.8, 0.3]
# Temp_P1_T1_29 = [462, 20.8, 9.15, 0.3]
# Temp_P1_T1_30 = [472, 21.9, 8.48, 0.3]
# Temp_P1_T1_31 = [484, 23.2, 6.40, 0.3]

data_2d = np.asarray([
    [76.0, 1.629, 20.75, 10], [82.5, 1.790, 21.30, 10], [95.5, 1.768, 20.50, 10],
    [104.8, 1.838, 20.66, 10], [110.6, 2.098, 21.93, 10], [116.5, 1.976, 23.86, 10],
    [125, 2.138, 21.21, 10], [138.5, 2.493, 20.25, 10], [145.5, 2.614, 22.93, 10],
    [155, 2.602, 20.59, 10], [163, 2.899, 20.59, 10], [170, 3.294, 20.71, 10],
    [175, 3.550, 20.15, 10], [180, 3.736, 19.21, 10], [186, 4.004, 18.51, 10],
    [200, 4.9, 16.1, 2.750], [210, 5.6, 19.2, 2.709], [220, 6.0, 19.0, 2.753],
    [229, 6.2, 15.5, 2.683], [237, 6.9, 16.8, 2.608], [245, 7.1, 18.9, 2.649],
    [252, 7.8, 18.2, 2.602], [260, 8.1, 16.0, 2.441], [266, 8.2, 16.5, 2.437],
    [274, 8.4, 16.2, 2.469], [284, 8.9, 17.0, 2.383], [294, 9.3, 14.2, 2.228],
    [305, 9.9, 16.6, 2.228], [318, 10.4, 17.6, 2.165], [331, 11.3, 15.1, 1.993],
    [353, 12.2, 12.2, 1.526], [364, 12.2, 14.9, 1.249], [375, 13.7, 13.3, 1.049],
    [386, 13.3, 12.2, 0.810], [396, 14.9, 13.5, 0.654], [404, 15.3, 13.4, 0.506],
    [412, 15.4, 12.8, 0.422], [424, 17.0, 11.7, 0.355], [438, 17.2, 11.4, 0.3],
    [449, 19.4, 10.8, 0.3], [462, 20.8, 9.15, 0.3], [472, 21.9, 8.48, 0.3],
    [484, 23.2, 6.40, 0.3]
])

# data_2d = np.asarray(
#     [Temp_P1_T1_1, Temp_P1_T1_2, Temp_P1_T1_3, Temp_P1_T1_4,
#      Temp_P1_T1_5, Temp_P1_T1_6, Temp_P1_T1_7, Temp_P1_T1_7,
#      Temp_P1_T1_8, Temp_P1_T1_9, Temp_P1_T1_10, Temp_P1_T1_10,
#      Temp_P1_T1_11, Temp_P1_T1_12, Temp_P1_T1_13, Temp_P1_T1_14,
#      Temp_P1_T1_15, Temp_P1_T1_16, Temp_P1_T1_17, Temp_P1_T1_18,
#      Temp_P1_T1_19, Temp_P1_T1_20, Temp_P1_T1_21, Temp_P1_T1_22,
#      Temp_P1_T1_23, Temp_P1_T1_24, Temp_P1_T1_25, Temp_P1_T1_26,
#      Temp_P1_T1_27, Temp_P1_T1_28, Temp_P1_T1_29, Temp_P1_T1_30,
#      Temp_P1_T1_31])

for i, d in enumerate(data_2d): # convert units to SI
    data_2d[i] = np.multiply(d, units)

# make data array
Temp_array = data_2d[:, 0]
P1_array = data_2d[:, 1]
T1_total_array = data_2d[:, 2]
PSD_array = data_2d[:, 3]

# extract QP induced part
T1_QP_array = [1/(1/T-1/T1_other) for T in T1_total_array]
Up_total_array = T1_total_array*(1/P1_array-1) #Upward transition rate

# convert time to rates
Gamma_down_total = [1/T for T in T1_total_array]
Gamma_down_QP = [1/T for T in T1_QP_array]
Gamma_up_QP = [1/T for T in Up_total_array]
Gamma_PSD = [1/T for T in PSD_array]
Gamma_pa = [getGamma_pa(Temp)+1/0.002750 for Temp in Temp_array]


Gamma_down_QP_scale = [np.log(1/G) for G in Gamma_down_QP]
Gamma_up_QP_scale = [np.log(1/G) for G in Gamma_up_QP]
Gamma_PSD_scale = [np.log(1/G) for G in Gamma_PSD]
Gamma_pa_scale = [np.log(1/G) for G in Gamma_pa]


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
# plt.plot(Temp_array, Gamma_up_QP, label='1/Up')
# plt.plot(Temp_array, Gamma_down_QP, label='1/T1_QP')
# plt.plot(Temp_array, Gamma_down_total, label='1/T1_tot')
# plt.plot(np.array([186, 305, 365, 425])*10**(-3), [766, 1701, 2184, 2891], label='QB Up Rate')
# plt.plot(Temp_array[15:38], Gamma_PSD[15:38], label='PSD Measured')
# plt.plot(Temp_array, Gamma_pa, label='Photon Assisted Prediction')
# plt.yscale("log")
# plt.xlabel('Temp (Kelvin)')
# plt.ylabel('Rate (Hz)')
# plt.grid()
# plt.legend()
# plt.show()

# """Plot Up/PSD ratio"""
# Up_PSD_ratio = np.divide(Gamma_up_QP, Gamma_PSD)
# plt.plot(Temp_array[:27], Up_PSD_ratio[:27], label='Up/PSD')
# plt.yscale("log")
# plt.xlabel('Temp (Kelvin)')
# plt.ylabel('Ratio (Up(P1, T1)/PSD)')
# plt.grid()
# plt.legend()
# plt.show()


"""Processed data and function fit"""
# T_Pre = np.array([186, 305, 365, 425])*10**(-3)
# Gamma_Pre = np.array([766, 1701, 2184, 2891])*1.0
T_Pre = np.array([236, 253, 267, 280, 293, 307, 318, 332, 343, 355, 366, 378,
                  403, 414, 426, 437, 451, 463, 475, 488, 501])\
        *10**(-3)
Gamma_Pre = np.array([874, 1065, 1251, 1325, 1191, 1391, 1518, 1723, 1681,
                      1696, 2029, 2051, 2369, 2970, 2878, 3444, 3411, 3568, 3915,
                      4061, 3952])\
            *1.0
# Gamma_pa_scale = [np.log(1/G) for G in Gamma_pa]

plt.scatter([1/T for T in T_Pre], [np.log(1/G) for G in Gamma_Pre], label='QB_Up_Pre')
# plt.scatter([1/T for T in T_Pre], [0.45/T-8.7 for T in T_Pre], label='Fit, f0=9.4GHz')
plt.plot([1/T for T in T_Pre], [0.45/T-8.75 for T in T_Pre], label='Fit, f0=9.4GHz')
# plt.scatter([1/T for T in T_Pre], [1.2/T-10.8 for T in T_Pre], label='Fit, f0=9.4GHz')
plt.plot([1/T for T in T_Pre], [1.2/T-10.8 for T in T_Pre], label='Fit, f0=25GHz')
# plt.scatter([1/T for T in Temp_array], Gamma_down_QP_scale, label='QP_down')
# plt.scatter([1/T for T in Temp_array], Gamma_up_QP_scale, label='QP_Up')
# plt.scatter([1/T for T in Temp_array], [0.06/T-7.5 for T in Temp_array], label='QP_fit')
# # plt.scatter([1/T for T in Temp_array], [0.5/T-10.5 for T in Temp_array], label='QP_fit')
# plt.scatter([1/T for T in Temp_array[:27]], Gamma_PSD_scale[:27], label='PSD')
# # plt.scatter([1/T for T in Temp_array], Gamma_pa_scale, label='Photon Assisted')
plt.xlabel('1/Temp(Kelvin)')
plt.ylabel('ln(df/Gamma)')
plt.grid()
plt.legend()
plt.show()