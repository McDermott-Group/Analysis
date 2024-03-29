### To analyze the P1 vs J2 Bias data
### src: Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorT1P1_2021Aug18_HighDensity\P1_J2_Q1

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorT1P1_2021Aug18_HighDensity'

file_Number = np.arange(0, 50, 1)

experiment_name_Q1 = ('P1_J2_Q1')
file_list_Q1 = [file_path.format(date, experiment_name_Q1, experiment_name_Q1) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_J2_Q1 = P1_JSweep()
P1_J2_Q1.add_data_from_matlab(file_list_Q1)

experiment_name_Q4 = ('P1_J2_Q4')
file_list_Q4 = [file_path.format(date, experiment_name_Q4, experiment_name_Q4) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_J2_Q4 = P1_JSweep()
P1_J2_Q4.add_data_from_matlab(file_list_Q4)

### J2 Bias offset
J2_Offset = 19.2
J2_Bias_Q1 = copy.deepcopy(P1_J2_Q1.J2_Bias)
J2_Bias_Q4 = copy.deepcopy(P1_J2_Q4.J2_Bias)
P1_J2_Q1 = copy.deepcopy(P1_J2_Q1.occ_1D_avg)
P1_J2_Q4 = copy.deepcopy(P1_J2_Q4.occ_1D_avg)
# print('J2_Bias=', J2_Bias)
# print('P1_J2_Q1=', P1_J2_Q1)
# print('P1_J2_Q4=', P1_J2_Q4)

J2_Bias_Q1 = 1000*J2_Bias_Q1-J2_Offset
J2_Bias_Q4 = 1000*J2_Bias_Q4-J2_Offset

print('J2_Bias_Q1=', J2_Bias_Q1)
print('J2_Bias_Q4=', J2_Bias_Q4)
print('P1_J2_Q1=', P1_J2_Q1)
print('P1_J2_Q4=', P1_J2_Q4)

# plt.plot(J2_Bias_Q1, P1_J2_Q1, label='Q1')
# plt.plot(J2_Bias_Q4, P1_J2_Q4, label='Q4')
# plt.xlabel('J2 Bias (mDAC)')
# plt.ylabel('P1')
# plt.grid()
# plt.legend()
# plt.show()

plt.plot(J2_Bias_Q1*4.8, P1_J2_Q1, label='Q1')
plt.plot(J2_Bias_Q4*4.8, P1_J2_Q4, label='Q4')
plt.xlabel('J2 Bias (GHz)')
plt.ylabel('P1')
plt.grid()
plt.legend()
plt.show()