"""
Used to analyze the up/down transition rate based on the T1, P1 as
a function of blackbody temperature
"""

import numpy as np
import matplotlib.pyplot as plt

# array for each temp in mK, P1 in units of %, T1 in units of
T1_other = 20.0
Temp_P1_T1_1 = [85, 2.75, 16.1]
Temp_P1_T1_2 = [127, 3, 16.3]
Temp_P1_T1_3 = [185, 4.5, 17.9]
Temp_P1_T1_4 = [238, 7.25, 14.7]
Temp_P1_T1_5 = [305, 11.25, 14.4]
Temp_P1_T1_6 = [365, 15, 14.1]
Temp_P1_T1_7 = [425, 18, 10.8]
Temp_P1_T1_8 = [485, 26, 6.6]
Temp_P1_T1_9 = [545, 33, 2.2]

data_2d = np.asarray([Temp_P1_T1_1, Temp_P1_T1_2, Temp_P1_T1_3, Temp_P1_T1_4,
           Temp_P1_T1_5, Temp_P1_T1_6, Temp_P1_T1_7, Temp_P1_T1_7,
           Temp_P1_T1_8, Temp_P1_T1_9])

Temp_array = data_2d[:, 0]
P1_array = data_2d[:, 1]
T1_total_array = data_2d[:, 2]
T1_QP_array = [1/(1/T-1/T1_other) for T in T1_total_array]
Tup_array = T1_total_array*(100/P1_array-1) #QP induced upward transition rate

# print('T1_total_array=', T1_total_array)
# print('T1_QP_array=', T1_QP_array)
# print('Tup_array=', Tup_array)

# plt.plot(Temp_array, P1_array)
# plt.xlabel('Temp (mK)')
plt.plot(Temp_array, Tup_array, label='Up')
plt.plot(Temp_array, T1_QP_array, label='T1_QP')
plt.plot(Temp_array, T1_total_array, label='T1_tot')
plt.yscale("log")
plt.xlabel('Temp (mK)')
plt.ylabel('Time (us)')
plt.grid()
plt.legend()
plt.show()