"""
PSD for Q1, Q2, Q4 at different poison power and recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\07-21-21
!!! The data set may be not the same thing since the repetition rate is different so the poisoning length is different
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q1_20us = np.array([[-90, 0.601], [-9, 0.559], [-6, 0.571], [-3, 0.577], [0, 0.443]])

Q2_10us = np.array([[-90, 45], [-10, 43.1], [0, 34.2], [3, 20.5], [6, 11.7]])
Q2_20us = np.array([[-90, 45], [-9, 45.9], [-6, 41.5], [-3, 41.9], [0, 34.7], [3, 21.0], [6, 13.0], [10, 1.62]])

Q4_10us = np.array([[-90, 2.92], [3, 1.85]])
Q4_20us = np.array([[-90, 2.92], [-9, 2.80], [-6, 2.89], [-3, 2.76], [0, 2.64], [3, 1.99], [6, 1.30]])

# Convert ms to Hz
Data_list = [Q1_20us, Q2_10us, Q2_20us, Q4_10us, Q4_20us]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]

plt.plot(Q1_20us[1:, 0], Q1_20us[1:, 1], label='Q1_20us')
plt.plot(Q2_10us[1:, 0], Q2_10us[1:, 1], label='Q2_10us')
plt.plot(Q2_20us[1:, 0], Q2_20us[1:, 1], label='Q2_20us')
plt.plot(Q4_10us[1:, 0], Q4_10us[1:, 1], label='Q4_10us')
plt.plot(Q4_20us[1:, 0], Q4_20us[1:, 1], label='Q4_20us')
plt.xlabel('Poison Power (dBm)')
plt.ylabel('QPT (Hz)')
# plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()
