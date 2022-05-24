### To analyze the P1 vs J2 Bias data
### src: Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorP1_2021Aug25_HighDensity\

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorP1_2021Aug26_HighDensity'

file_Number = np.arange(0, 200, 1)

### J2 Bias offset
experiment_name_Q2 = ('P1_J2_Q2')
file_list_Q2 = [file_path.format(date, experiment_name_Q2, experiment_name_Q2) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q2 = P1_JSweep()
P1_Q2.add_data_from_matlab(file_list_Q2)


### J2 Bias offset
# J2_Offset = 19.2
J2_Bias_Q2 = copy.deepcopy(P1_Q2.J2_Bias)
P1_J2_Q2 = copy.deepcopy(P1_Q2.occ_1D_avg)
Std_J2_Q2 = copy.deepcopy(P1_Q2.occ_1D_std)
J2_Bias_Q2 = 1000*J2_Bias_Q2
print('J2_Bias_Q2=', J2_Bias_Q2)
print('P1_J2_Q2=', P1_J2_Q2)
print('std_J2_Q2=', Std_J2_Q2)

plt.errorbar(J2_Bias_Q2, P1_J2_Q2, yerr=Std_J2_Q2, label='Q2')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

