"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')

date = '2021Mar17PSD'
power = 'Neg6'

experiment_name_PSD = ('Interleave_PSD_'+power)
PSD_file_Number = np.arange(0, 10, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]


QPT_Q4 = QPTunneling_Liu(name='Q4 with poison={}'.format(power))
QPT_Q4.add_datasets(PSD_file, HMM=True)

QPT_List = [QPT_Q4]
plotMultiFittedPSD(QPT_List)

p1_psd = [QPT_Q4.T_parity]
print(power, p1_psd)
