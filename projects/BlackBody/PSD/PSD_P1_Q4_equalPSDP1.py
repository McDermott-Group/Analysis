"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = '04-15-21'
# temp = 'Neg80'
# date = '04-16-21'
# temp = '83mK'
# temp = 'Neg305'
date = 'PSDP1T1Study2021Apr17'
temp = '266mK'
P1 = 0.5

experiment_name_P1 = ('Interleave_P1_'+temp)
experiment_name_PSD = ('Interleave_PSD_'+temp)
P1_file_Number = np.arange(0, 5, 1)
PSD_file_Number = P1_file_Number
P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
P1CleanDirty = OneStateCleanDirty()
P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=P1)
PSD_Clean_Files = P1CleanDirty.clean_files

QPT_Q4_Poison_Clean = QPTunneling_Liu(name='P1={:.3f} at temp={}'.format(P1CleanDirty.clean_P1, temp))
# QPT_Q4_Poison_Clean = QPTunneling_Wilen(name='P1={:.3f} at temp={}'.format(P1CleanDirty.clean_P1, temp))
QPT_Q4_Poison_Clean.add_datasets(PSD_Clean_Files)

QPT_List = [QPT_Q4_Poison_Clean]
plotMultiFittedPSD(QPT_List)

p1_psd = [P1CleanDirty.clean_P1, QPT_Q4_Poison_Clean.T_parity]
print (temp, p1_psd)