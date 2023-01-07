### To analyze the P1 vs J2 Bias data
### src: Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorP1_2021Aug25_HighDensity\

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = '09-28-21'

# file_Number = np.arange(43, 53, 1)
# file_Number = np.arange(53, 63, 1)
# file_Number = np.arange(5, 40, 1)
file_Number = np.arange(0, 10, 1)

### J2 Bias offset
experiment_name_Q1 = ('P1_JB_Q2')
file_list_Q1 = [file_path.format(date, experiment_name_Q1, experiment_name_Q1) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q1 = P1_JSweep()
P1_Q1.add_data_from_matlab(file_list_Q1)


### J2 Bias offset
# J2_Offset = 19.2
J2_Bias_Q1 = copy.deepcopy(P1_Q1.J2_Bias)
P1_J2_Q1 = copy.deepcopy(P1_Q1.occ_1D_avg)
Std_J2_Q1 = copy.deepcopy(P1_Q1.occ_1D_std)
J2_Bias_Q1 = 1000*J2_Bias_Q1
# print('J2_Bias_Q1=', J2_Bias_Q1)
# print('P1_J2_Q1=', P1_J2_Q1)
# print('std_J2_Q1=', Std_J2_Q1)

# plt.errorbar(J2_Bias_Q1, P1_J2_Q1, yerr=Std_J2_Q1, label='Q1')
plt.errorbar(J2_Bias_Q1, P1_J2_Q1, label='Q2')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

