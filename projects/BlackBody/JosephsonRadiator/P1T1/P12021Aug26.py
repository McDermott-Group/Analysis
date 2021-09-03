### To analyze the P1 vs J2 Bias data
### src: Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorP1_2021Aug25_HighDensity\

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorP1_2021Aug25_HighDensity'

file_Number = np.arange(0, 50, 1)

experiment_name_Q1 = ('P1_J2_Q1')
file_list_Q1 = [file_path.format(date, experiment_name_Q1, experiment_name_Q1) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q1 = P1_JSweep()
P1_Q1.add_data_from_matlab(file_list_Q1)

experiment_name_Q2 = ('P1_J2_Q2')
file_list_Q2 = [file_path.format(date, experiment_name_Q2, experiment_name_Q2) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q2 = P1_JSweep()
P1_Q2.add_data_from_matlab(file_list_Q2)

experiment_name_Q3 = ('P1_J2_Q3')
file_list_Q3 = [file_path.format(date, experiment_name_Q3, experiment_name_Q3) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q3 = P1_JSweep()
P1_Q3.add_data_from_matlab(file_list_Q3)

experiment_name_Q4 = ('P1_J2_Q4')
file_list_Q4 = [file_path.format(date, experiment_name_Q4, experiment_name_Q4) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q4 = P1_JSweep()
P1_Q4.add_data_from_matlab(file_list_Q4)

### J2 Bias offset
# J2_Offset = 19.2
J2_Bias_Q1 = copy.deepcopy(P1_Q1.J2_Bias)
J2_Bias_Q2 = copy.deepcopy(P1_Q2.J2_Bias)
J2_Bias_Q3 = copy.deepcopy(P1_Q3.J2_Bias)
J2_Bias_Q4 = copy.deepcopy(P1_Q4.J2_Bias)
P1_J2_Q1 = copy.deepcopy(P1_Q1.occ_1D_avg)
P1_J2_Q2 = copy.deepcopy(P1_Q2.occ_1D_avg)
P1_J2_Q3 = copy.deepcopy(P1_Q3.occ_1D_avg)
P1_J2_Q4 = copy.deepcopy(P1_Q4.occ_1D_avg)
Std_J2_Q1 = copy.deepcopy(P1_Q1.occ_1D_std)
Std_J2_Q2 = copy.deepcopy(P1_Q2.occ_1D_std)
Std_J2_Q3 = copy.deepcopy(P1_Q3.occ_1D_std)
Std_J2_Q4 = copy.deepcopy(P1_Q4.occ_1D_std)
# print('J2_Bias=', J2_Bias)
# print('P1_J2_Q1=', P1_J2_Q1)
# print('P1_J2_Q4=', P1_J2_Q4)

J2_Bias_Q1 = 1000*J2_Bias_Q1
J2_Bias_Q2 = 1000*J2_Bias_Q2
J2_Bias_Q3 = 1000*J2_Bias_Q3
J2_Bias_Q4 = 1000*J2_Bias_Q4

print('J2_Bias_Q1=', J2_Bias_Q1)
print('J2_Bias_Q2=', J2_Bias_Q2)
print('J2_Bias_Q3=', J2_Bias_Q3)
print('J2_Bias_Q4=', J2_Bias_Q4)
print('P1_J2_Q1=', P1_J2_Q1)
print('P1_J2_Q2=', P1_J2_Q2)
print('P1_J2_Q3=', P1_J2_Q3)
print('P1_J2_Q4=', P1_J2_Q4)
print('std_J2_Q1=', Std_J2_Q1)
print('std_J2_Q2=', Std_J2_Q2)
print('std_J2_Q3=', Std_J2_Q3)
print('std_J2_Q4=', Std_J2_Q4)

plt.errorbar(J2_Bias_Q1, P1_J2_Q1, yerr=Std_J2_Q1, label='Q1')
plt.errorbar(J2_Bias_Q2, P1_J2_Q2, yerr=Std_J2_Q2, label='Q2')
plt.errorbar(J2_Bias_Q3, P1_J2_Q3, yerr=Std_J2_Q3, label='Q3')
plt.errorbar(J2_Bias_Q4, P1_J2_Q4, yerr=Std_J2_Q4, label='Q4')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

# plt.plot(J2_Bias_Q1*4.8, P1_J2_Q1, label='Q1')
# plt.plot(J2_Bias_Q4*4.8, P1_J2_Q4, label='Q2')
# plt.xlabel('J2 Bias (GHz)')
# plt.ylabel('P1')
# plt.grid()
# plt.legend()
# plt.show()