"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')

date = '02-17-21'
power = 'C400mD_'
P1 = 0.1

experiment_name_P1 = (power+'Interleave_P1_Neg100')
experiment_name_PSD = (power+'Interleave_PSD_Neg100')
P1_file_Number = np.arange(60, 90, 1)
PSD_file_Number = P1_file_Number
P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
P1CleanDirty = OneStateCleanDirty()
P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=P1)
PSD_Clean_Files = P1CleanDirty.clean_files
PSD_Dirty_Files = P1CleanDirty.dirty_files

QPT_Q4_Poison_Clean = QPTunneling_Liu(name='Q4_Clean_P1<{} with poison={}'.format(P1, power))
QPT_Q4_Poison_Dirty = QPTunneling_Liu(name='Q4_Dirty_P1>{} with poison={}'.format(P1, power))
QPT_Q4_Poison_Clean.add_datasets(PSD_Clean_Files)
QPT_Q4_Poison_Dirty.add_datasets(PSD_Dirty_Files)

# experiment_name_P1 = ('Interleave_P1_Neg10')
# experiment_name_PSD = ('Interleave_PSD_Neg10')
#
# P1_file_Number = np.arange(0, 100, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.04)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg10dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg10_Clean (P1<0.04)'+date)
# QPT_Q4_Poison_Neg10dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg10_Dirty (P1>0.04)'+date)
#
# QPT_Q4_Poison_Neg10dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg10dBm_Dirty.add_datasets(PSD_Dirty_Files)

### New data after Thanks giving fridge cycling

QPT_List = [QPT_Q4_Poison_Clean, QPT_Q4_Poison_Dirty]
#
plotMultiFittedPSD(QPT_List)
