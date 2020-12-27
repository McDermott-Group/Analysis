"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')

# date = '10-29-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)

# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Off_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Off_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

date = '10-30-20'
# experiment_name_P1 = ('Interleave_P1_Neg20')
# experiment_name_PSD = ('Interleave_PSD_Neg20')
#
# P1_file_Number = np.arange(0, 250, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.04)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg20dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Clean (P1<0.04)')
# QPT_Q4_Poison_Neg20dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Dirty (P1>0.04)')

# QPT_Q4_Poison_Neg20dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg20dBm_Dirty.add_datasets(PSD_Dirty_Files)


# experiment_name_P1 = ('Interleave_P1_Neg16')
# experiment_name_PSD = ('Interleave_PSD_Neg16')
#
# P1_file_Number = np.arange(0, 250, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.04)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg16dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Clean (P1<0.04)')
# QPT_Q4_Poison_Neg16dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Dirty (P1>0.04)')
#
# QPT_Q4_Poison_Neg16dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg16dBm_Dirty.add_datasets(PSD_Dirty_Files)

# experiment_name_P1 = ('Interleave_P1_Neg13')
# experiment_name_PSD = ('Interleave_PSD_Neg13')
#
# P1_file_Number = np.arange(0, 250, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.05)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg13dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg13dBm_Clean (P1<0.05)')
# QPT_Q4_Poison_Neg13dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg13dBm_Dirty (P1>0.05)')
#
# QPT_Q4_Poison_Neg13dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg13dBm_Dirty.add_datasets(PSD_Dirty_Files)

# experiment_name_P1 = ('Interleave_P1_Neg10')
# experiment_name_PSD = ('Interleave_PSD_Neg10')
#
# P1_file_Number = np.arange(0, 1000, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.075)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg10dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg10dBm_Clean (P1<0.075)')
# QPT_Q4_Poison_Neg10dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg10dBm_Dirty (P1>0.075)')
#
# QPT_Q4_Poison_Neg10dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg10dBm_Dirty.add_datasets(PSD_Dirty_Files)

"""Q6"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')

# date = '12-07-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.07)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Off_Clean (P1<0.07)')
# QPT_Q6_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Off_Dirty (P1>0.07)')
#
# QPT_Q6_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-07-20'
# experiment_name_P1 = ('Interleave_P1_Neg10')
# experiment_name_PSD = ('Interleave_PSD_Neg10')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.09)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg10dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg10dBm_Clean (P1<0.09)')
# QPT_Q6_Poison_Neg10dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg10dBm_Dirty (P1>0.09)')
#
# QPT_Q6_Poison_Neg10dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg10dBm_Dirty.add_datasets(PSD_Dirty_Files)


# date = '12-08-20'
# experiment_name_P1 = ('Interleave_P1_Neg20')
# experiment_name_PSD = ('Interleave_PSD_Neg20')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.08)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg20dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg20dBm_Clean (P1<0.08)')
# QPT_Q6_Poison_Neg20dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg20dBm_Dirty (P1>0.08)')
#
# QPT_Q6_Poison_Neg20dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg20dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-09-20'
# experiment_name_P1 = ('Interleave_P1_Neg16')
# experiment_name_PSD = ('Interleave_PSD_Neg16')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.08)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg16dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg16dBm_Clean (P1<0.08)')
# QPT_Q6_Poison_Neg16dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg16dBm_Dirty (P1>0.08)')
#
# QPT_Q6_Poison_Neg16dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg16dBm_Dirty.add_datasets(PSD_Dirty_Files)

date = '12-09-20'
experiment_name_P1 = ('Interleave_P1_Neg13')
experiment_name_PSD = ('Interleave_PSD_Neg13')

P1_file_Number = np.arange(0, 500, 1)
PSD_file_Number = P1_file_Number
P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
P1CleanDirty = OneStateCleanDirty()
P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.085)

PSD_Clean_Files = P1CleanDirty.clean_files
PSD_Dirty_Files = P1CleanDirty.dirty_files

QPT_Q6_Poison_Neg13dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg13dBm_Clean (P1<0.085)')
QPT_Q6_Poison_Neg13dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg13dBm_Dirty (P1>0.085)')

QPT_Q6_Poison_Neg13dBm_Clean.add_datasets(PSD_Clean_Files)
QPT_Q6_Poison_Neg13dBm_Dirty.add_datasets(PSD_Dirty_Files)

# QPT_List = [QPT_Q4_Poison_Neg100dBm_Clean, QPT_Q4_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg20dBm_Clean, QPT_Q4_Poison_Neg20dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg16dBm_Clean, QPT_Q4_Poison_Neg16dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg13dBm_Clean, QPT_Q4_Poison_Neg13dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg10dBm_Clean, QPT_Q4_Poison_Neg10dBm_Dirty]

# QPT_List = [QPT_Q6_Poison_Neg100dBm_Clean, QPT_Q6_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg20dBm_Clean, QPT_Q6_Poison_Neg20dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg16dBm_Clean, QPT_Q6_Poison_Neg16dBm_Dirty]
QPT_List = [QPT_Q6_Poison_Neg13dBm_Clean, QPT_Q6_Poison_Neg13dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg10dBm_Clean, QPT_Q6_Poison_Neg10dBm_Dirty]

plotMultiFittedPSD(QPT_List)
