### To analyze the P1 vs Duration
### Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\10-27-21

from antennalib import P1_Avg_vs_Any
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/Antenna/SUXmon2/LIU/VitoChip2/{}/{}/MATLABData/{}')

date = '10-28-21'
file_Number = np.arange(0, 5, 1)

QB_Name = 'Q4'
experiment_name = ('{}_J7_P1_FixedPoison'.format(QB_Name))
file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

P1 = P1_Avg_vs_Any()
P1.add_data_from_matlab(file_list, data_type2='J7_Idle')

P1_2D = []

for i in range(len(P1.var)):
    var = P1.var[i]
    p1 = int(P1.p1_1D_avg[i]*100000)/100000.0
    std = int(P1.p1_1D_std[i]*100000)/100000.0
    P1_2D.append([var, p1, std])

# J_P1_2D = np.array(J_P1_2D)
print('QB=', QB_Name)
print('P1_2D=', P1_2D)

plt.plot(P1.var, P1.p1_1D_avg, label=QB_Name)
plt.xlabel(P1.var_name+' (us)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()