"""
Poison at Q1, measure PSD for Q2 at same poison power and length but different recover time
Raw data directory
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\Q1Poison_Q4Recover_2021Jul27

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import copy
from scipy.constants import *
# from antennalib import getGamma_pa

# = [[Power1, Rate1], [Power2, Rate2]] Power is the RF output power (dBm), rate is ms
Q4_3dBm = np.array([
    [10, 0.93042], [20, 1.09279],  [30, 1.36891],  [40, 1.38267], [50, 1.5126],
    [60, 1.75090], [70, 1.70607], [80, 1.79082], [90, 1.82861], [100, 1.97047],
    [110, 2.01873], [120, 2.01715], [130, 2.09584], [140, 2.17907], [150, 2.24128],
    [160, 2.22086], [170, 2.28237], [180, 2.30563], [190, 2.37257], [200, 2.44782],
    [210, 2.49027], [220, 2.49559]])
# Convert ms to Hz
Data_list = [Q4_3dBm]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]

### for fitting to get time constant


def f_exp_off(t_recover, amp, T1, off):
    return amp * np.exp(-t_recover / T1) + off

time = Q4_3dBm[:, 0]
PSD = Q4_3dBm[:, 1]
params, covariance = curve_fit(f_exp_off, time, PSD)
amp, T1, off = params[0], params[1], params[2]
# print(amp, T1, off)

### plot
plt.plot(Q4_3dBm[:, 0], Q4_3dBm[:, 1], label='Poisoning 3dB higher than 268uV')
plt.plot(Q4_3dBm[:, 0], [amp * np.exp(-t / T1) + off for t in time],
         '--', label='fit time constant = {} (us)'.format(T1))

plt.xlabel('Recover after poison (us)')
plt.ylabel('Q2 QPT (Hz)')
# plt.xscale('log')
plt.grid()
plt.legend()
plt.title('Q4 Parity Rate vs Recover Time After Q1 Poison')
plt.show()

