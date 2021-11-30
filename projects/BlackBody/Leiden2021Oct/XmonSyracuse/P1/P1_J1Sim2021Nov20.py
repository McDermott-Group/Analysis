### To analyze the P1 vs JJ Bias data
### Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\2021Nov20_P1_J1Sim


from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')
date = '2021Nov20_P1_J1Sim'

file_Number_all = np.arange(0, 89, 1)
f_highP1 = [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 28, 29, 30, 34, 35, 36, 37,
            38, 47, 48, 49, 51, 65, 66, 67, 70, 71, 72, 80, 81, 82, 83, 85, 86, 87]

# file_Number_all = np.concatenate((file_Number_all_1, file_Number_all_1))
# f_highP1 = np.concatenate((f_highP1_1, f_highP1_2))

file_Number = [x for x in file_Number_all if x not in f_highP1]

J_P1_2D = []

QB_Name = 'Q1'

for k in range(6):
    experiment_name = ('P1_J1Slow_{}_Chunk{}'.format(QB_Name, str(k)))
    file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

    # experiment_name = ('P1_J1Slow_{}'.format(QB_Name))
    # file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

    P1 = P1_JSweep()
    P1.add_data_from_matlab(file_list, data_type2='J1_Slow_Bias')
    # J_Bias_P1 = copy.deepcopy(P1.J2_Bias)

    for i in range(len(P1.J_Bias)):
        J_Bias = int(P1.J_Bias[i]*10000)/10000.0
        p1 = int(P1.occ_1D_avg[i]*100000)/100000.0
        std = int(P1.occ_1D_std[i]*100000)/100000.0
        J_P1_2D.append([J_Bias, p1, std])

J_P1_2D = np.array(J_P1_2D)
print('QB=', QB_Name)
print('J_P1_2D=', J_P1_2D)

# plt.plot(P1.J_Bias, P1.occ_1D_avg, label=QB_Name)
# # plt.plot(P1.J_Bias, P1.occ_1D_avg[::-1], label='Reversed')
# plt.xlabel('J1 Slow Bias (mV)')
# plt.ylabel('P1')
# # plt.yscale('log')
# plt.grid()
# plt.legend()
# plt.show()