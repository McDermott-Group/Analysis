"""
2021Oct
Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev01\Leiden_2020Feb\LIU\Q1\03-17-20\QP_Tunneling_PSD\MATLABData
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev01/Leiden_2020Feb/LIU/Q1/{}/{}/MATLABData/{}')

date = '03-16-20'
experiment_name_PSD = ('QP_Tunneling_PSD_2')
PSD_file_Number = np.arange(1, 30, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
QPT_Q = QPTunneling_Wilen(name='Gap')
# QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
#                           format(QB_id, str(JB), str(JB*4.6)))
QPT_Q.add_datasets(PSD_file)
QPT_List = [QPT_Q]
plotMultiFittedPSD(QPT_List, save=False, name='Gap')

QPT_Q.get_fit()
# print(QB_id)
# print(J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
rate = QPT_Q.params[0]  # units is Hz
rate = int(rate*100)/100.0
# print('rate=', rate)
J_QPT_2D.append([JB, rate])
# print('J_QPT_2D=', J_QPT_2D)

print('QB_id=', QB_id)
print('J_QPT_2D=', J_QPT_2D)


