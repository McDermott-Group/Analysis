"""
Poison at Q1, measure PSD for Q2 at same poison power and length but different recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\Q1Poison_Q2Recover_2021Jul26

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import copy
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q2_3dBm = np.array([
    [10, 1.74665], [20, 2.30689], [30, 2.9389], [40, 3.33774], [50, 3.94282],
    [60, 4.67817], [70, 4.68675], [80, 5.25995], [90, 5.29739], [100, 5.91467],
    [110, 6.19221], [120, 6.38775], [130, 7.13684], [140, 7.39163], [150, 7.75057],
    [160, 8.25281], [170, 8.30428], [180, 8.98687], [190, 8.8382], [200, 9.29845],
    [210, 9.81104], [220, 9.67608]])

# Convert ms to Hz
Data_list = [Q2_3dBm]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]

### for fitting to get time constant


def f_exp_off(t_recover, amp, T1, off):
    return amp * np.exp(-t_recover / T1) + off

time = Q2_3dBm[:, 0]
PSD = Q2_3dBm[:, 1]
params, covariance = curve_fit(f_exp_off, time, PSD)
amp, T1, off = params[0], params[1], params[2]
# print(amp, T1, off)

### plot
plt.plot(Q2_3dBm[:, 0], Q2_3dBm[:, 1], label='Poisoning 3dB higher than 268uV')
plt.plot(Q2_3dBm[:, 0], [amp * np.exp(-t / T1) + off for t in time],
         '--', label='fit time constant = {} (us)'.format(T1))

plt.xlabel('Recover after poison (us)')
plt.ylabel('Q2 QPT (Hz)')
plt.xscale('log')
plt.grid()
plt.legend()
plt.title('Q2 Parity Rate vs Recover Time After Q1 Poison')
plt.show()

