### To analyze the P1 vs JJ Bias data
### src: Z:\mcdermott-group\data\Antenna\Circmon\Liu\CW20180514A_Ox2\10-15-21\P1_J6_Q1
### Z:\mcdermott-group\data\Antenna\Circmon\Liu\CW20180514A_Ox2\2021Oct16_QB124_P1_J6Radiator
### Z:\mcdermott-group\data\Antenna_Temporary\Circmon\Liu\CW20180514A_Ox2\2021Nov02_QB124_P1_J6Radiator_Sim


from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
# date = '2021Oct16_QB124_P1_J6Radiator'
# date = '2021Nov02_QB124_P1_J6Radiator_Sim'
date = '2021Nov04_QB124_P1PSD_J8Radiator_Sim'
file_Number = np.arange(0, 50, 1)

J_P1_2D = []

QB_Name = 'Q2'

for k in range(5):
    print('k=', k)
    experiment_name = ('P1_J8Slow_{}_Chunk{}'.format(QB_Name, str(k)))
    file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1 = P1_JSweep()
    P1.add_data_from_matlab(file_list, data_type2='J8_Slow_Bias')
    # J_Bias_P1 = copy.deepcopy(P1.J2_Bias)

    for i in range(len(P1.J_Bias)):
        J_Bias = int(P1.J_Bias[i]*10000)/10000.0
        p1 = int(P1.occ_1D_avg[i]*100000)/100000.0
        std = int(P1.occ_1D_std[i]*100000)/100000.0
        J_P1_2D.append([J_Bias, p1, std])

# J_P1_2D = np.array(J_P1_2D)
print('QB=', QB_Name)
print('J_P1_2D=', J_P1_2D)

# plt.plot(J_P1_2D[:, 0], J_P1_2D[:, 1], label=QB_Name)
# plt.xlabel('J6 Bias (mDAC)')
# plt.ylabel('P1')
# # plt.yscale('log')
# plt.grid()
# plt.legend()
# plt.show()