"""
Used to analyze the up/down transition rate based on the T1, P1 as
a function of blackbody temperature
"""

from antennalib import getGamma_pa

import numpy as np
import matplotlib.pyplot as plt

"""Processed data and function fit"""
T_Pre_Q4 = np.array([236, 253, 267, 280, 293, 307, 318, 332, 343, 355, 366, 378,
                  403, 414, 426, 437, 451, 463, 475, 488, 501])\
        *10**(-3)
Gamma_Pre_Q4 = np.array([874, 1065, 1251, 1325, 1191, 1391, 1518, 1723, 1681,
                      1696, 2029, 2051, 2369, 2970, 2878, 3444, 3411, 3568, 3915,
                      4061, 3952])\
            *1.0

T_Pre = np.array([78, 185, 237, 274, 306, 336, 365, 397, 426])\
        *10**(-3)

Gamma_Pre_Q1 = np.array([359, 334, 903, 1670, 1990, 1559, 2455, 2637, 1679])\
            *1.0

Gamma_Pre_Q2 = np.array([560, 174, 590, 842, 960, 1203, 1292, 1422, 1777])\
            *1.0

Gamma_Pre_Q3 = np.array([1664, 1538, 2641, 2798, 3147, 2894, 3496, 3481, 4018])\
            *1.0


# plt.scatter([1/T for T in T_Pre], [np.log(1/G) for G in Gamma_Pre_Q1], label='Q1_Up_Pre')
# plt.scatter([1/T for T in T_Pre], [np.log(1/G) for G in Gamma_Pre_Q2], label='Q2_Up_Pre')
plt.scatter([1/T for T in T_Pre], [np.log(1/G) for G in Gamma_Pre_Q3], label='Q3_Up_Pre')
# plt.plot([1/T for T in T_Pre_Q4], [np.log(1/G) for G in Gamma_Pre_Q4], label='Q4_Up_Pre')
plt.plot([1/T for T in T_Pre], [0.2/T-8.7 for T in T_Pre], label='Fit, f0=4GHz')
plt.xlabel('1/Temp(Kelvin)')
plt.ylabel('ln(df/Gamma)')
plt.grid()
plt.legend()
plt.show()