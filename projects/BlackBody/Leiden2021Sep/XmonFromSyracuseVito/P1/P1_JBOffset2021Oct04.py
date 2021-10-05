### To analyze the P1 vs J7 Bias data
### src: Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\10-04-21\P1_JB_Q3\MATLABData\

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')
date = '10-04-21'
file_Number = np.arange(2, 12, 1)

experiment_name_Q1 = ('P1_JB_Q3')
file_list_Q1 = [file_path.format(date, experiment_name_Q1, experiment_name_Q1) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_Q1 = P1_JSweep()
P1_Q1.add_data_from_matlab(file_list_Q1, data_type2='JB7_Bias')
J2_Bias_Q1 = copy.deepcopy(P1_Q1.J2_Bias)

plt.plot(J2_Bias_Q1, P1_Q1.occ_1D_avg, label='Q1')
plt.xlabel('J7 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()