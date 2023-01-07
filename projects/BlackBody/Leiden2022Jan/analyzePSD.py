"""
2021Oct04
### src: Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\10-04-21\
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

QP_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')

date = '05-16-22'
QB_id = 'Q3'

PSD_file_Number = np.arange(5, 15, 1)
# PSD_file_Number = np.arange(45, 70, 1)
# PSD_file_Number = np.arange(30, 40, 1)

experiment_name_PSD = (QB_id+'_PSD')
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
QPT_Q = QPTunneling_Wilen(name='{}'.
                          format(QB_id))
QPT_Q.add_datasets(PSD_file)
QPT_List = [QPT_Q]
plotMultiFittedPSD(QPT_List, save=False, name='{}'.
                          format(QB_id))

QPT_Q.get_fit()
rate = QPT_Q.params[0]  # units is Hz
rate = int(rate*100)/100.0
# print('rate=', rate)
# J_QPT_2D.append([rate])
# print('J_QPT_2D=', J_QPT_2D)

# print('QB_id=', QB_id)
# print('J_QPT_2D=', J_QPT_2D)