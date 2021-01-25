"""
2020Dec03
This is for clean and dirty state
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
#
# date = '10-29-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Off_Clean (P1<0.035) Oct29')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Off_Dirty (P1>0.035) Oct29')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '10-30-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
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

# date = '12-09-20'
# experiment_name_P1 = ('Interleave_P1_Neg13')
# experiment_name_PSD = ('Interleave_PSD_Neg13')
#
# P1_file_Number = np.arange(0, 500, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.085)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q6_Poison_Neg13dBm_Clean = QPTunneling_Wilen(name='Q6_Poison_Neg13dBm_Clean (P1<0.085)')
# QPT_Q6_Poison_Neg13dBm_Dirty = QPTunneling_Wilen(name='Q6_Poison_Neg13dBm_Dirty (P1>0.085)')
#
# QPT_Q6_Poison_Neg13dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q6_Poison_Neg13dBm_Dirty.add_datasets(PSD_Dirty_Files)

### New data after Thanks giving fridge cycling

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
# date = '12-23-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.030)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Clean (P1<0.030) Dec23')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Dirty (P1>0.030) Dec23')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-28-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.030)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Clean (P1<0.030) Dec28')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Dirty (P1>0.030) Dec28')
# #
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-30-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.030)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Clean (P1<0.030) Dec30')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Dirty (P1>0.030) Dec30')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-30-20'
# experiment_name_P1 = ('Interleave_P1_Neg20')
# experiment_name_PSD = ('Interleave_PSD_Neg20')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg20dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg20dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg20dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg20dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg19')
# experiment_name_PSD = ('Interleave_PSD_Neg19')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 281, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg19dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg19dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg19dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg19dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg19dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg19dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg18')
# experiment_name_PSD = ('Interleave_PSD_Neg18')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg18dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg18dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg18dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg18dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg18dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg18dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg17')
# experiment_name_PSD = ('Interleave_PSD_Neg17')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg17dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg17dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg17dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg17dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg17dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg17dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg16')
# experiment_name_PSD = ('Interleave_PSD_Neg16')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg16dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg16dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg16dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg16dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg15')
# experiment_name_PSD = ('Interleave_PSD_Neg15')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg15dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg15dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg15dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg15dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg15dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg15dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg14')
# experiment_name_PSD = ('Interleave_PSD_Neg14')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.035)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg14dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg14dBm_Clean (P1<0.035)')
# QPT_Q4_Poison_Neg14dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg14dBm_Dirty (P1>0.035)')
#
# QPT_Q4_Poison_Neg14dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg14dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-31-20'
# experiment_name_P1 = ('Interleave_P1_Neg13')
# experiment_name_PSD = ('Interleave_PSD_Neg13')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 243, 1)
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

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg12')
# experiment_name_PSD = ('Interleave_PSD_Neg12')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.05)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg12dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg12dBm_Clean (P1<0.05)')
# QPT_Q4_Poison_Neg12dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg12dBm_Dirty (P1>0.05)')
#
# QPT_Q4_Poison_Neg12dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg12dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg11')
# experiment_name_PSD = ('Interleave_PSD_Neg11')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.05)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg11dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg11dBm_Clean (P1<0.05)')
# QPT_Q4_Poison_Neg11dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg11dBm_Dirty (P1>0.05)')
#
# QPT_Q4_Poison_Neg11dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg11dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg10')
# experiment_name_PSD = ('Interleave_PSD_Neg10')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(100, 400, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.060)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg10dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg10dBm_Clean (P1<0.060)')
# QPT_Q4_Poison_Neg10dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg10dBm_Dirty (P1>0.060)')
#
# QPT_Q4_Poison_Neg10dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg10dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg9')
# experiment_name_PSD = ('Interleave_PSD_Neg9')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.070)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg9dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg9dBm_Clean (P1<0.070)')
# QPT_Q4_Poison_Neg9dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg9dBm_Dirty (P1>0.070)')
#
# QPT_Q4_Poison_Neg9dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg9dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg8')
# experiment_name_PSD = ('Interleave_PSD_Neg8')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.080)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg8dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg8dBm_Clean (P1<0.080)')
# QPT_Q4_Poison_Neg8dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg8dBm_Dirty (P1>0.080)')
#
# QPT_Q4_Poison_Neg8dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg8dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg7')
# experiment_name_PSD = ('Interleave_PSD_Neg7')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.10)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg7dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg7dBm_Clean (P1<0.10)')
# QPT_Q4_Poison_Neg7dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg7dBm_Dirty (P1>0.10)')
#
# QPT_Q4_Poison_Neg7dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg7dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-01-21'
# experiment_name_P1 = ('Interleave_P1_Neg6')
# experiment_name_PSD = ('Interleave_PSD_Neg6')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 300, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.10)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg6dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg6dBm_Clean (P1<0.10)')
# QPT_Q4_Poison_Neg6dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg6dBm_Dirty (P1>0.10)')
#
# QPT_Q4_Poison_Neg6dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg6dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-07-21'
# experiment_name_P1 = ('Remove_Interleave_P1_Neg100')
# experiment_name_PSD = ('Remove_Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 10, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.029)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Clean (P1<0.030) Dec23')
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q4_Poison_Neg100dBm_Dirty (P1>0.030) Dec23')
#
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

"""Q2"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q2_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-06-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.10)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q2_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q2_Poison_Neg100dBm_Clean (P1<0.10) Nov06')
# QPT_Q2_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q2_Poison_Neg100dBm_Dirty (P1>0.10) Nov06')
#
# QPT_Q2_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q2_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '01-06-21'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.10)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q2_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q2_Poison_Neg100dBm_Clean (P1<0.10) Jan06')
# QPT_Q2_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q2_Poison_Neg100dBm_Dirty (P1>0.10) Jan06')
#
# QPT_Q2_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q2_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

"""Q1"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q1_withQ2Poison/{}/{}/MATLABData/{}')
# date = '11-16-20'
# experiment_name_P1 = ('Interleave_P1_Neg100')
# experiment_name_PSD = ('Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.1)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q1_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q1_Poison_Neg100dBm_Clean (P1<0.1) Nov16')
# QPT_Q1_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q1_Poison_Neg100dBm_Dirty (P1>0.1) Nov16')
#
# QPT_Q1_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q1_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

# date = '12-10-20'
# experiment_name_P1 = ('Charge_Interleave_P1_Neg100')
# experiment_name_PSD = ('Charge_Interleave_PSD_Neg100')
# print('power=', experiment_name_P1)
#
# P1_file_Number = np.arange(0, 200, 1)
# PSD_file_Number = P1_file_Number
# P1_file = [QP_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]
# P1CleanDirty = OneStateCleanDirty()
# P1CleanDirty.update_experiment_name(experiment_name_P1, experiment_name_PSD)
# P1CleanDirty.add_p1_data_from_matlab(P1_file, P1_threshold=0.10)
#
# PSD_Clean_Files = P1CleanDirty.clean_files
# PSD_Dirty_Files = P1CleanDirty.dirty_files
#
# QPT_Q1_Poison_Neg100dBm_Clean = QPTunneling_Wilen(name='Q1_Poison_Neg100dBm_Clean (P1<0.10) Dec10')
# QPT_Q1_Poison_Neg100dBm_Dirty = QPTunneling_Wilen(name='Q1_Poison_Neg100dBm_Dirty (P1>0.10) Dec10')
#
# QPT_Q1_Poison_Neg100dBm_Clean.add_datasets(PSD_Clean_Files)
# QPT_Q1_Poison_Neg100dBm_Dirty.add_datasets(PSD_Dirty_Files)

QPT_List = [QPT_Q4_Poison_Neg100dBm_Clean, QPT_Q4_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg20dBm_Clean, QPT_Q4_Poison_Neg20dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg19dBm_Clean, QPT_Q4_Poison_Neg19dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg18dBm_Clean, QPT_Q4_Poison_Neg18dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg17dBm_Clean, QPT_Q4_Poison_Neg17dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg16dBm_Clean, QPT_Q4_Poison_Neg16dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg15dBm_Clean, QPT_Q4_Poison_Neg15dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg14dBm_Clean, QPT_Q4_Poison_Neg14dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg13dBm_Clean, QPT_Q4_Poison_Neg13dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg12dBm_Clean, QPT_Q4_Poison_Neg12dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg11dBm_Clean, QPT_Q4_Poison_Neg11dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg10dBm_Clean, QPT_Q4_Poison_Neg10dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg9dBm_Clean, QPT_Q4_Poison_Neg9dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg8dBm_Clean, QPT_Q4_Poison_Neg8dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg7dBm_Clean, QPT_Q4_Poison_Neg7dBm_Dirty]
# QPT_List = [QPT_Q4_Poison_Neg6dBm_Clean, QPT_Q4_Poison_Neg6dBm_Dirty]

# QPT_List = [QPT_Q6_Poison_Neg100dBm_Clean, QPT_Q6_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg20dBm_Clean, QPT_Q6_Poison_Neg20dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg16dBm_Clean, QPT_Q6_Poison_Neg16dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg13dBm_Clean, QPT_Q6_Poison_Neg13dBm_Dirty]
# QPT_List = [QPT_Q6_Poison_Neg10dBm_Clean, QPT_Q6_Poison_Neg10dBm_Dirty]

# QPT_List = [QPT_Q2_Poison_Neg100dBm_Clean, QPT_Q2_Poison_Neg100dBm_Dirty]

# QPT_List = [QPT_Q1_Poison_Neg100dBm_Clean, QPT_Q1_Poison_Neg100dBm_Dirty]
#
plotMultiFittedPSD(QPT_List)
