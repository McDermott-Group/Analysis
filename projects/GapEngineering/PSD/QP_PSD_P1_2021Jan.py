"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

# """Q1"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q1_withQ2Poison/{}/{}/MATLABData/{}')
#
# date = '02-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 30, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.01)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q1_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q1_Poison_Off_Clean (P1<0.01)')
# QPT_Q1_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q1_Poison_Off_Dirty (P1>0.01)')
#
# QPT_Q1_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q1_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

"""Q2"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q2_withQ5Poison/{}/{}/MATLABData/{}')
#
# date = '02-02-21'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 40, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.055)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q2_Poison_Neg100dBm_Clean = QPTunneling_Liu(name='Q2_Poison_Off_Clean (P1<0.055)')
# QPT_Q2_Poison_Neg100dBm_Dirty = QPTunneling_Liu(name='Q2_Poison_Off_Dirty (P1>0.055)')
#
# QPT_Q2_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q2_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)


# """Q3"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q3_withQ2Poison/{}/{}/MATLABData/{}')
#
# date = '02-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.008)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q3_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q3_Poison_Off_Clean (P1<0.008)')
# QPT_Q3_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q3_Poison_Off_Dirty (P1>0.008)')
#
# QPT_Q3_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q3_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)


"""Q4"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
#
# date = '02-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 20, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.02)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Off_Clean (P1<0.02)')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Off_Dirty (P1>0.02)')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

"""Q6"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2021Jan/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')

# date = '01-28-21'
# date = '02-03-21'
date = '2021Feb06'
# date = '02-08-21'
# date = '2021Feb09'
experiment_name_P1 = ('Interleave_P1_Neg100')
experiment_name_PSD = ('Interleave_PSD_Neg100')

P1_file_Number = np.arange(0, 100, 1)
PSD_file_Number = P1_file_Number
P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
P1CleanDirty = OneStateCleanDirty()
P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.02)

PSD_Clean_Files = P1CleanDirty.clean_files
PSD_Dirty_Files = P1CleanDirty.dirty_files

QPT_Q6_Poison_Neg100dBm_Clean = QPTunneling_Liu(name='Q6_Poison_Off_Clean (P1<0.02)'+date)
QPT_Q6_Poison_Neg100dBm_Dirty = QPTunneling_Liu(name='Q6_Poison_Off_Dirty (P1>0.02)'+date)

QPT_Q6_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
QPT_Q6_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '02-05-21'
# experiment_name_P1 = ('Interleave_P1_Neg10')
# experiment_name_PSD = ('Interleave_PSD_Neg10')
#
# P1_file_Number = np.arange(0, 100, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.1)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg10dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg10dBm_Clean (P1<0.1)')
# QPT_Q6_Poison_Neg10dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg10dBm_Dirty (P1>0.1)')
#
# QPT_Q6_Poison_Neg10dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg10dBm_Dirty.add_datasets(PSD_Dirty_Files)

### New data after Thanks giving fridge cycling

# QPT_List = [QPT_Q1_Poison_Neg100dBm_Clean, QPT_Q1_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q2_Poison_Neg100dBm_Clean, QPT_Q2_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q3_Poison_Neg100dBm_Clean, QPT_Q3_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg100dBm_Clean, QPT_Q4_Poison_Neg100dBm_Dirty]
QPT_List = [QPT_Q6_Poison_Neg100dBm_Clean, QPT_Q6_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg20dBm_Clean, QPT_Q6_Poison_Neg20dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg16dBm_Clean, QPT_Q6_Poison_Neg16dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg13dBm_Clean, QPT_Q6_Poison_Neg13dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg10dBm_Clean, QPT_Q6_Poison_Neg10dBm_Dirty]

# QPT_List = [QPT_Q2_Poison_Neg100dBm_Clean, QPT_Q2_Poison_Neg100dBm_Dirty]

# QPT_List = [QPT_Q1_Poison_Neg100dBm_Clean, QPT_Q1_Poison_Neg100dBm_Dirty]
#
plotMultiFittedPSD(QPT_List)
