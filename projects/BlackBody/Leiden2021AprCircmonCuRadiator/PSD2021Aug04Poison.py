"""
PSD for Q1, Q2, Q4 at different poison power and recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\PSD_Poison_R3Q3_2021Aug04

"""
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q2_20us = np.array(
    [[-90, 42.96027], [-6, 41.66562], [-3, 41.21706], [0, 41.39949], [1, 36.56876],
     [2, 30.82847], [3, 28.11687], [4, 23.38326], [5, 22.83866], [6, 20.28172], [7, 16.98884]])

Q4_20us = np.array(
    [[-90, 2.75153], [-6, 2.81873], [-3, 2.73644], [0, 2.71434], [1, 2.49889],
     [2, 2.25679], [3, 2.10055], [4, 1.64581], [5, 1.39076], [6, 1.24397], [7, 1.09749]]
)

# Convert ms to Hz
Data_list = [Q2_20us, Q4_20us]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]
# print('Q2_20us=', Q2_20us)



### plot
plt.plot(Q2_20us[1:, 0], Q2_20us[1:, 1], label='Q2_20us')
plt.plot(Q4_20us[1:, 0], Q4_20us[1:, 1], label='Q4_20us')

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
