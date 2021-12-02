"""
PSD for Q1, Q2, Q4 at different J6 Bias
Data:
Z:\mcdermott-group\data\Antenna\Circmon\Liu\CW20180514A_Ox2\2021Oct16_QB124_PSD_J6Radiator
Fitting method Chris' no white noise verison

"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np

# [f_QB(GHz), 2df(MHz), Ej(GHz), QP_Lifetime(ms)]
Q1 = np.array([
    [4.602, 1.4, 26*0.345, 1.16], [4.510, 1.8, 25.12*0.345, 0.833], [4.371, 2.5, 23.6*0.345, 0.55]
])

Q4 = np.array([
    [4.941, 1.1, 27*0.363, 5.6], [4.638, 1.9, 24.5*0.363, 4.03], [4.217, 7, 20*0.363, 3.16]
])

Q1[:, -1] = 1000/Q1[:, -1]
Q4[:, -1] = 1000/Q4[:, -1]


plt.plot(Q1[:, -2], Q1[:, -1], color='b', label='Q1')
plt.plot(Q4[:, -2], Q4[:, -1], color='y', label='Q4')
# plt.xlabel('QB Charge Disperison (MHz)')
plt.xlabel('EJ(GHz)')
plt.ylabel('PSD (Hz)')
plt.xscale('log')
plt.yscale('log')
plt.grid(True, which="both")
plt.legend(loc=1)
plt.show()




