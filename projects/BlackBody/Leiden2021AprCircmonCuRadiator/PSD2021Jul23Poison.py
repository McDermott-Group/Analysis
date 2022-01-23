"""
PSD for Q2, Q4 at different poison power and recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\PSD_Poison_R1Q1_2021Jul23

"""
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q2_20us = np.array([[-90, 45], [-9, 38.15], [-6, 30.5], [-3, 12.99], [0, 6.94], [3, 2.36], [6, 0.837]])

Q4_20us = np.array([[-90, 2.92], [-9, 2.778], [-6, 2.687], [-3, 2.334], [0, 1.674], [3, 1.033], [6, 0.550], [9, 0.328]])

# Convert ms to Hz
Data_list = [Q2_20us, Q4_20us]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]
# print('Q2_20us=', Q2_20us)



### plot
# plt.plot(Q2_20us[1:, 0], Q2_20us[1:, 1], label='Q2_20us')
# plt.plot(Q4_20us[1:, 0], Q4_20us[1:, 1], label='Q4_20us')

# plt.hlines(Q2_20us[0][1], Q2_20us[1][0], Q2_20us[-1][0], linestyles='dashed', label='Q2_Base')
# plt.hlines(Q4_20us[0][1], Q4_20us[1][0], Q4_20us[-1][0], linestyles='dashed', label='Q4_Base')

plt.hlines(Q2_20us[0][1], Q2_20us[1][0], Q2_20us[-1][0], linestyles='dashed', label='No Poisoning Baseline')
plt.plot(Q2_20us[1:, 0], Q2_20us[1:, 1], label='With Poisoning')

plt.xlabel('Q1 Poison Power (dBm)')
plt.ylabel('Q2 QPT (Hz)')
# plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()
plt.title('Q2 Parity Rate vs Q1 Poisoning Power')
plt.show()
# plt.draw()

### the added part from Q3
# Extract the added channel
# Q2_20us_added = copy.deepcopy(Q2_20us)
# Q4_20us_added = copy.deepcopy(Q4_20us)
#
# Data_Added_list = [Q2_20us_added, Q4_20us_added]
# for d in Data_Added_list:
#     d[:, 1] = d[:, 1] - d[0][1]
#
# # print('Q2_20us_added=', Q2_20us_added)
#
# plt.figure()
# plt.plot(Q2_20us_added[1:, 0], Q2_20us_added[1:, 1], label='Q2_20us_added')
# plt.plot(Q4_20us_added[1:, 0], Q4_20us_added[1:, 1], label='Q4_20us_added')
# plt.xlabel('Poison Power (dBm)')
# plt.ylabel('QPT (Hz)')
# plt.yscale('log')
# plt.grid()
# plt.legend()
# plt.show()
