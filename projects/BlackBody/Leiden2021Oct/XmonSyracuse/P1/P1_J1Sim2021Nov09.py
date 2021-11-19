### To analyze the P1 vs JJ Bias data
### Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\10-20-21\P1_J7_Q2\MATLABData
### Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\P12021Oct22_Q1WithJ7
### Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\10-27-21\P1_J1_Q2\MATLABData
### Z:\mcdermott-group\data\Antenna_Temporary\SUXmon2\Liu\VitoChip2\11-02-21
### Z:\mcdermott-group\data\Antenna_Temporary\Circmon\Liu\CW20180514A_Ox2\11-02-21\P1_J6Slow_Q4
### Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\2021Nov08_Q3_PSD_J1Radiator\P1_J1Slow_Q3\MATLABData

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')

# date = '11-04-21'
date = '2021Nov08_Q3_PSD_J1Radiator'
# date = 'P12021Oct26_Q123WithJ7'
file_Number = np.arange(0, 50, 1)

J_P1_2D = []

QB_Name = 'Q3'

experiment_name = ('P1_J1Slow_{}'.format(QB_Name))
file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

P1 = P1_JSweep()
P1.add_data_from_matlab(file_list, data_type2='J1_Slow_Bias')

for i in range(len(P1.J_Bias)):
    J_Bias = int(P1.J_Bias[i]*10000)/10000.0
    p1 = int(P1.occ_1D_avg[i]*100000)/100000.0
    std = int(P1.occ_1D_std[i]*100000)/100000.0
    J_P1_2D.append([J_Bias, p1, std])

# J_P1_2D = np.array(J_P1_2D)
print('QB=', QB_Name)
print('J_P1_2D=', J_P1_2D)

plt.plot(P1.J_Bias, P1.occ_1D_avg, label=QB_Name)
# plt.plot(P1.J_Bias, P1.occ_1D_avg[::-1]-0.005, label='Reversed')
# plt.xlabel('J1 Bias (mDAC)')
plt.xlabel('J1 Bias (mV)')

# plt.plot(4604*P1.J_Bias, P1.occ_1D_avg, label=QB_Name)
# plt.xlabel('J7 Bias (GHz)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()