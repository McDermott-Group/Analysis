"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q6"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')

date = '2021Feb24'
power = 'Neg17'

experiment_name_PSD = ('Long_Interleave_PSD_'+power)
PSD_file_Number = np.arange(0, 2, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in P1_file_Number]


QPT_Q6_Poison_Clean = QPTunneling_Wilen(name='Q6_Clean_P1<{} with poison={}'.format(P1, power))
QPT_Q6_Poison_Clean.add_datasets(PSD_file)

QPT_List = [QPT_Q6_Poison_Clean]
plotMultiFittedPSD(QPT_List)

p1_psd = [P1CleanDirty.clean_P1, QPT_Q6_Poison_Clean.T_parity]
print (power, p1_psd)