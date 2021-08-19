"""
PSD for Q1, Q2, Q4 at different poison power and recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\PSD_Poison_R3Q3_2021Jul22

"""
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q1_20us = np.array([[-90, 0.59258], [-9, 0.50993], [-6, 0.49965], [-3, 0.31055], [0, 0.42209], [3, 0.26633]])

Q2_20us = np.array([[-90, 44.276], [-9, 44.75776], [-6, 44.07476], [-3, 44.0434], [0, 39.18505], [3, 29.144], [6, 20.1006], [9, 8.37631]])

Q4_20us = np.array([[-90, ]])

# Convert ms to Hz
Data_list = [Q1_20us, Q2_20us, Q4_20us]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]
# print('Q2_20us=', Q2_20us)



### plot
plt.plot(Q1_20us[1:, 0], Q1_20us[1:, 1], label='Q1_20us')
plt.plot(Q2_20us[1:, 0], Q2_20us[1:, 1], label='Q2_20us')
plt.plot(Q4_20us[1:, 0], Q4_20us[1:, 1], label='Q4_20us')

plt.hlines(Q1_20us[0][1], Q1_20us[1][0], Q1_20us[-1][0], linestyles='dashed', label='Q1_Base')
plt.hlines(Q2_20us[0][1], Q2_20us[1][0], Q2_20us[-1][0], linestyles='dashed', label='Q2_Base')
plt.hlines(Q4_20us[0][1], Q4_20us[1][0], Q4_20us[-1][0], linestyles='dashed', label='Q4_Base')

plt.xlabel('Poison Power (dBm)')
plt.ylabel('QPT (Hz)')
# plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()
# plt.show()
plt.draw()

### the added part from Q3
# Extract the added channel
Q2_20us_added = copy.deepcopy(Q2_20us)
Q4_20us_added = copy.deepcopy(Q4_20us)

Data_Added_list = [Q2_20us_added, Q4_20us_added]
for d in Data_Added_list:
    d[:, 1] = d[:, 1] - d[0][1]

# print('Q2_20us_added=', Q2_20us_added)

plt.figure()
plt.plot(Q2_20us_added[1:, 0], Q2_20us_added[1:, 1], label='Q2_20us_added')
plt.plot(Q4_20us_added[1:, 0], Q4_20us_added[1:, 1], label='Q4_20us_added')
plt.xlabel('Poison Power (dBm)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()
