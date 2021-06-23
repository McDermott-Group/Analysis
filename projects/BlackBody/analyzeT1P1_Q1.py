"""
Used to analyze the up/down transition rate based on the T1, P1 as
a function of blackbody temperature. This data is Q4
"""

from antennalib import getGamma_pa

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *

Temp_array = [72, 79, 94, 103, 116, 124, 138, 145, 150, 160, 170, 178, 186,
              193, 200, 205, 210, 216, 220, 225, 229, 234, 238]
Temp_array = [t*1.0/1000 for t in Temp_array]
Gamma_up = np.array([1677.642002132063, 1457.3107713146997, 1385.2256294069566,
            1533.7428603145759,
            1826.368116440605, 1908.936285507801, 1970.287911592732,
            1996.6553207550817,
            1782.1123859683423, 2220.591857277547, 2519.35636046541,
            2527.0192248417834,
            2676.148910089545, 3474.349367624638, 3084.695547065121,
            3449.771927464184,
            3496.763095270301, 3704.311882306011, 3646.2012835151704,
            3467.4735325469746,
            3377.0245594049393, 3549.1613414160825, 3775.9947643572787])

Gamma_up_Error =np.array([92.08281661311533, 90.16374952793146, 66.96708030704261,
                           49.6545535620637, 130.33465459476025, 150.97276051858975,
                           110.3285186997869, 105.45651006581053, 28.49962120556832,
                           105.57729533704818, 115.50174430919787, 77.26132977135258,
                           116.03838248715147, 216.960799774514, 82.57574612317244,
                           89.73950870898395, 118.12520292808043, 165.92415035722385,
                           113.47157958961809, 49.15345359566372, 85.38184726680979,
                           103.56974165372954, 88.71666834810173])
GammaErrorMinus = Gamma_up - Gamma_up_Error
GammaErrorPlus = Gamma_up + Gamma_up_Error
Gamma_up_fit = Temp_array
Gamma_up_fit_2 = Temp_array
Gamma_up_fit_3 = Temp_array
Temp_kelvin = Temp_array
# for i in range(len(Temp_array)):
#     Temp_kelvin[i] = Temp_array[i]*1.0/1000
#     df1 = 5e3
#     df2 = 5e3
#     df3 = 5e3
#     f1 = 5
#     f2 = 6
#     f3 = 6
#     Gamma_up_fit[i] = df1 / (np.exp(f1 / (20.8 * Temp_kelvin[i]) - 1))
    # Gamma_up_fit_2[i] = df2 / (np.exp(f2 / (20.8 * Temp_kelvin[i]) - 1))
    # Gamma_up_fit_3[i] = df3 / (np.exp(f3 / (20.8 * Temp_kelvin[i]) - 1))

# print(Temp_kelvin)

"""Plot rate"""
# plt.plot(Temp_array, Gamma_up, label='1/Up')
# plt.plot(Temp_array, Gamma_up_fit, label='1/Up_fit')
# plt.plot(Temp_array, Gamma_up_fit_2, label='1/Up_fit_2')
# plt.plot(Temp_array, Gamma_up_fit_3, label='1/Up_fit_3')
# plt.xlabel('Temp (Kelvin)')
# plt.ylabel('Rate (Hz)')
# plt.grid()
# plt.legend()
# plt.show()

# """Plot Up/PSD ratio"""
plt.scatter([1/T for T in Temp_array], [np.log(1/G) for G in Gamma_up], label='QB_Up_Pre')
plt.fill_between([1 / T for T in Temp_array],
                 [np.log(1 / G) for G in GammaErrorMinus],
                 [np.log(1 / G) for G in GammaErrorPlus], alpha=0.2)
plt.plot([1/T for T in Temp_array], [0.03/T-7.7 for T in Temp_array], label='Fit, f0=0.624GHz')
plt.plot([1/T for T in Temp_array], [0.22/T-9.1 for T in Temp_array], label='Fit, f0=4.57GHz')
# plt.scatter([1/T for T in T_Pre], [1.2/T-10.8 for T in T_Pre], label='Fit, f0=9.4GHz')
plt.xlabel('1/Temp(Kelvin)')
plt.ylabel('ln(df/Gamma)')
plt.grid()
plt.legend()
plt.show()
