### To analyze the P1 vs JJ Bias data
### Z:\mcdermott-group\data\Antenna_Temporary\SUXmon2\Liu\VitoChip2\P12021Oct31_Q2WithJ7OnChip_TimeDomain\P1_J7_Q2_Chunk0\MATLABData

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna_Temporary/SUXmon2/LIU/VitoChip2/{}/{}/MATLABData/{}')
date = 'P12021Oct31_Q2WithJ7OnChip_TimeDomain'
file_Number = np.arange(0, 50, 1)

J_P1_2D = []

QB_Name = 'Q4'

for k in range(5):
    experiment_name = ('P1_J7_{}_Chunk{}'.format(QB_Name, str(k)))
    file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

# experiment_name = ('P1_J7_{}'.format(QB_Name))
# file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

    P1 = P1_JSweep()
    P1.add_data_from_matlab(file_list, data_type2='J7_Bias')
    # J_Bias_P1 = copy.deepcopy(P1.J2_Bias)

    for i in range(len(P1.J_Bias)):
        J_Bias = int(P1.J_Bias[i]*10000)/10000.0
        p1 = int(P1.occ_1D_avg[i]*100000)/100000.0
        std = int(P1.occ_1D_std[i]*100000)/100000.0
        J_P1_2D.append([J_Bias, p1, std])

# J_P1_2D = np.array(J_P1_2D)
print('QB=', QB_Name)
print('J_P1_2D=', J_P1_2D)

# plt.plot(P1.J_Bias, P1.occ_1D_avg, label=QB_Name)
plt.plot(J_P1_2D[:, 0], J_P1_2D[:, 0], label=QB_Name)
# plt.plot(P1.J_Bias, P1.occ_1D_avg[::-1], label='Reversed')
plt.xlabel('J7 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()