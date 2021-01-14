from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np


QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
date = '10-29-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file_1 = np.arange(0, 25, 1)     # 473 Hz
# QP_file_1 = np.arange(30, 35, 1)    # 386 Hz
# QP_file_1 = np.arange(36, 45, 1)    # 499 Hz
# QP_file_1 = np.arange(70, 80, 1)    # 497 Hz
# QP_file_1 = np.arange(80, 84, 1)    # 357 Hz
# QP_file_1 = np.arange(90, 102, 1)    # 514 Hz
# QP_file_1 = np.arange(106, 120, 1)    # 500 Hz
# QP_file_1 = np.arange(150, 180, 1)    # 488 Hz
# QP_file_1 = np.arange(207, 210, 1)    # 400 Hz
# QP_file_1 = np.arange(226, 250, 1)    # 459 Hz
# QP_file_1 = np.arange(252, 255, 1)    # 375 Hz

# date = '10-29-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file_1 = np.arange(100, 120, 1)     # 384 Hz
# QP_file_1 = np.arange(165, 170, 1)     # 358 Hz
# QP_file_1 = np.arange(0, 50, 1)     # 390 Hz

# date = '10-29-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file_Clean = np.arange(10, 28, 1)
# QP_file_Clean = np.append(QP_file_Clean, np.arange(32, 42, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(98, 108, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(180, 200, 1))
#
# QP_file_Dirty = np.arange(29, 30, 1)
# # QP_file_Dirty = np.append(QP_file_Dirty, np.arange(44, 46, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(53, 59, 1))
# # QP_file_Dirty = np.append(QP_file_Dirty, np.arange(110, 111, 1))
# # QP_file_Dirty = np.append(QP_file_Dirty, np.arange(168, 169, 1))
# # QP_file_Dirty = np.append(QP_file_Dirty, np.arange(176, 177, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(427, 429, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(482, 490, 1))
#
# print('len(clean)=', len(QP_file_Clean))
# print('len(dirty)=', len(QP_file_Dirty))
#
# QP_filenames_Q4_Poison_Neg100dBm_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q4_Poison_Neg100dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg100dBm_Clean')
# QPT_Q4_Poison_Neg100dBm_Clean.add_datasets(QP_filenames_Q4_Poison_Neg100dBm_Clean, HMM=False)

#
# QP_filenames_Q4_Poison_Neg100dBm_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q4_Poison_Neg100dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg100dBm_Dirty')
# QPT_Q4_Poison_Neg100dBm_Dirty.add_datasets(QP_filenames_Q4_Poison_Neg100dBm_Dirty, HMM=False)

date = '10-30-20'
experiment_name = ('Interleave_PSD_Neg10')
# QP_file_1 = np.arange(0, 200, 1)     #  380 Hz
# QP_file_1 = np.arange(0, 8, 1)     #  341 Hz, WA 243 Hz
# QP_file_1 = np.arange(8, 59, 1)     # 390 Hz, WA 319 Hz
# QP_file_1 = np.arange(10, 50, 1)     # 400 Hz good
# QP_file_1 = np.arange(59, 64, 1)     # 353 Hz good
# QP_file_1 = np.arange(59, 60, 1)     # 428 Hz
# QP_file_1 = np.arange(60, 61, 1)     # 318 Hz
# QP_file_1 = np.arange(61, 62, 1)     # Bad fit and data 987 Hz
# QP_file_1 = np.arange(62, 63, 1)     # Bad fit and data 989 Hz
# QP_file_1 = np.arange(75, 95, 1)     # 464 Hz good
# QP_file_1 = np.arange(120, 125, 1)     # 329 Hz good
# QP_file_1 = np.arange(150, 175, 1)     # 411 Hz good
# QP_file_1 = np.arange(180, 185, 1)     # 269 Hz good
# QP_file_1 = np.arange(190, 206, 1)     # 464 Hz good
# QP_file_1 = np.arange(209, 217, 1)     # 370 Hz good
# QP_file_1 = np.arange(541, 555, 1)     # 370 Hz good
# QP_filenames_Q4_Poison_Neg10dBm = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_1]
# QPT_Q4_Poison_Neg10dBm = QPTunneling_Liu(name='Q4_Poison_Neg10dBm')
# QPT_Q4_Poison_Neg10dBm.add_datasets(QP_filenames_Q4_Poison_Neg10dBm, HMM=False)

# date = '10-30-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file_Clean = np.arange(0, 7, 1)
# QP_file_Clean = np.append(QP_file_Clean, np.arange(60, 63, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(180, 184, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(210, 215, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(282, 287, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(480, 486, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(502, 507, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(540, 555, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(560, 566, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(580, 586, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(630, 636, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(650, 656, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(700, 706, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(782, 786, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(890, 897, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(980, 985, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(1050, 1057, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(1091, 1097, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(1119, 1126, 1))
#
# QP_file_Dirty = np.arange(10, 55, 1)
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(75, 95, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(130, 145, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(150, 175, 1))
#
# print('len(clean)=', len(QP_file_Clean))
# print('len(dirty)=', len(QP_file_Dirty))
#
# QP_filenames_Q4_Poison_Neg10dBm_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q4_Poison_Neg10dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg10dBm_Clean')
# QPT_Q4_Poison_Neg10dBm_Clean.add_datasets(QP_filenames_Q4_Poison_Neg10dBm_Clean, HMM=False)
#
# QP_filenames_Q4_Poison_Neg10dBm_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q4_Poison_Neg10dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg10dBm_Dirty')
# QPT_Q4_Poison_Neg10dBm_Dirty.add_datasets(QP_filenames_Q4_Poison_Neg10dBm_Dirty, HMM=False)

# date = '10-30-20'
# experiment_name = ('Interleave_PSD_Neg13')
#
# QP_file_Clean = np.arange(20, 29, 1)
# QP_file_Clean = np.append(QP_file_Clean, np.arange(43, 49, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(102, 108, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(166, 169, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(220, 229, 1))
#
# QP_file_Dirty = np.arange(0, 7, 1)
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(10, 19, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(30, 38, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(60, 67, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(70, 80, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(90, 100, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(110, 118, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(130, 160, 1))
#
# print('len(clean)=', len(QP_file_Clean))
# print('len(dirty)=', len(QP_file_Dirty))
#
# QP_filenames_Q4_Poison_Neg13dBm_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q4_Poison_Neg13dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg13dBm_Clean')
# QPT_Q4_Poison_Neg13dBm_Clean.add_datasets(QP_filenames_Q4_Poison_Neg13dBm_Clean, HMM=False)
#
# QP_filenames_Q4_Poison_Neg13dBm_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q4_Poison_Neg13dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg13dBm_Dirty')
# QPT_Q4_Poison_Neg13dBm_Dirty.add_datasets(QP_filenames_Q4_Poison_Neg13dBm_Dirty, HMM=False)

# date = '10-30-20'
# experiment_name = ('Interleave_PSD_Neg16')
# QP_file_Clean = np.arange(27, 29, 1)
# QP_file_Clean = np.append(QP_file_Clean, np.arange(41, 48, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(67, 69, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(112, 119, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(164, 170, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(194, 200, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(233, 239, 1))
#
# QP_file_Dirty = np.arange(20, 26, 1)
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(30, 40, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(60, 66, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(90, 100, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(120, 130, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(150, 160, 1))
#
# print('len(clean)=', len(QP_file_Clean))
# print('len(dirty)=', len(QP_file_Dirty))
#
# QP_filenames_Q4_Poison_Neg16dBm_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q4_Poison_Neg16dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Clean')
# QPT_Q4_Poison_Neg16dBm_Clean.add_datasets(QP_filenames_Q4_Poison_Neg16dBm_Clean, HMM=False)
#
# QP_filenames_Q4_Poison_Neg16dBm_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q4_Poison_Neg16dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg16dBm_Dirty')
# QPT_Q4_Poison_Neg16dBm_Dirty.add_datasets(QP_filenames_Q4_Poison_Neg16dBm_Dirty, HMM=False)

#
# date = '10-30-20'
# experiment_name = ('Interleave_PSD_Neg20')
# QP_file_Clean = np.arange(10, 20, 1)
# QP_file_Clean = np.append(QP_file_Clean, np.arange(32, 50, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(90, 100, 1))
# QP_file_Clean = np.append(QP_file_Clean, np.arange(106, 113, 1))
#
# QP_file_Dirty = np.arange(5, 9, 1)
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(102, 105, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(134, 135, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(145, 149, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(183, 188, 1))
# QP_file_Dirty = np.append(QP_file_Dirty, np.arange(212, 216, 1))
#
# print('len(clean)=', len(QP_file_Clean))
# print('len(dirty)=', len(QP_file_Dirty))
#
# QP_filenames_Q4_Poison_Neg20dBm_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q4_Poison_Neg20dBm_Clean = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Clean')
# QPT_Q4_Poison_Neg20dBm_Clean.add_datasets(QP_filenames_Q4_Poison_Neg20dBm_Clean, HMM=False)
#
# QP_filenames_Q4_Poison_Neg20dBm_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q4_Poison_Neg20dBm_Dirty = QPTunneling_Liu(name='Q4_Poison_Neg20dBm_Dirty')
# QPT_Q4_Poison_Neg20dBm_Dirty.add_datasets(QP_filenames_Q4_Poison_Neg20dBm_Dirty, HMM=False)

date = '10-28-20'
experiment_name = ('Interleave_PSD_Neg100')
QP_file = np.arange(0, 100, 1)

QP_filenames_Q4_Poison_Neg100dBm = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
QPT_Q4_Poison_Neg100dBm = QPTunneling_Liu(name='Q4_Poison_Neg100dBm')
QPT_Q4_Poison_Neg100dBm.add_datasets(QP_filenames_Q4_Poison_Neg100dBm, HMM=False)

"""Q6 Starts"""

QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-02-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file_1 = np.arange(0, 10, 1)     # 298 Hz
# QP_file_1 = np.arange(0, 100, 1)     # 282 Hz

# QP_file_Clean = np.arange(65, 70, 1)     # 198 Hz Good
# QP_file_Clean = np.arange(442, 449, 1)     # 208 Good
# QP_file_Clean = np.arange(160, 170, 1)     # 234 Good
# QP_filenames_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q6_NoPoison_Clean = QPTunneling(name='Q6_NoPoisonClean')
# QPT_Q6_NoPoison_Clean.add_datasets(QP_filenames_Clean)

# QP_file_Dirty = np.arange(25, 30, 1)     # falied
# QP_file_Dirty = np.arange(70, 75, 1)     # 274 Hz
# QP_file_Dirty = np.arange(70, 90, 1)     # 290 Hz good
# QP_file_Dirty = np.arange(180, 200, 1)     # 306 Hz good
# QP_file_Dirty = np.arange(340, 350, 1)     #  failed
# QP_filenames_Dirty = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Dirty]
# QPT_Q6_NoPoison_Dirty = QPTunneling(name='Q6_NoPoisonDirty')
# QPT_Q6_NoPoison_Dirty.add_datasets(QP_filenames_Dirty)
# QPT.plot_psd()

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-03-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 100, 1)     # 221 Hz
# # QP_file = np.arange(0, 100, 1)     # 225 Hz
# # QP_file = np.arange(100, 150, 1)     # 234 Hz
# # QP_file = np.arange(210, 220, 1)     # 238 Hz
#
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_NoPoison = QPTunneling_Liu(name='Q6_NoPoison')
# QPT_Q6_NoPoison.add_datasets(QP_filenames)
#
# date = '11-03-20'
# experiment_name = ('Interleave_PSD_Neg20')
# QP_file = np.arange(0, 100, 1)     #
# QP_filenames_Neg20 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg20 = QPTunneling_Liu(name='Q6_Neg20')
# QPT_Q6_Neg20.add_datasets(QP_filenames_Neg20)
#
# date = '11-03-20'
# experiment_name = ('Interleave_PSD_Neg16')
# QP_file = np.arange(0, 100, 1)     #
# QP_filenames_Neg16 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg16 = QPTunneling_Liu(name='Q6_Neg16')
# QPT_Q6_Neg16.add_datasets(QP_filenames_Neg16)
#
# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg13')
# QP_file = np.arange(0, 100, 1)     #
# QP_filenames_Neg13 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg13 = QPTunneling_Liu(name='Q6_Neg13')
# QPT_Q6_Neg13.add_datasets(QP_filenames_Neg13)
#
# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 100, 1)     #
# QP_filenames_Neg10 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames_Neg10)


# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 400, 1)     # 221 Hz
# # QP_file = np.arange(0, 100, 1)     # 225 Hz
# # QP_file = np.arange(100, 150, 1)     # 234 Hz
# # QP_file = np.arange(210, 220, 1)     # 238 Hz
#
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_NoPoison = QPTunneling_Liu(name='Q6_NoPoison')
# QPT_Q6_NoPoison.add_datasets(QP_filenames)
#
# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg20')
# QP_file = np.arange(0, 400, 1)     #
# QP_filenames_Neg20 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg20 = QPTunneling_Liu(name='Q6_Neg20')
# QPT_Q6_Neg20.add_datasets(QP_filenames_Neg20)
#
# date = '11-05-20'
# experiment_name = ('Interleave_PSD_Neg16')
# QP_file = np.arange(0, 400, 1)     #
# QP_filenames_Neg16 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg16 = QPTunneling_Liu(name='Q6_Neg16')
# QPT_Q6_Neg16.add_datasets(QP_filenames_Neg16)
#
# date = '11-05-20'
# experiment_name = ('Interleave_PSD_Neg13')
# QP_file = np.arange(0, 400, 1)     #
# QP_filenames_Neg13 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg13 = QPTunneling_Liu(name='Q6_Neg13')
# QPT_Q6_Neg13.add_datasets(QP_filenames_Neg13)
#
# date = '11-05-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 400, 1)     #
# QP_filenames_Neg10 = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames_Neg10)


# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg100')
# # QP_file = np.arange(0, 1, 1)     #
# # QP_file = np.arange(0, 10, 1)     # 234 Hz
# QP_file = np.arange(34, 40, 1)     # 329 Hz (HMM 308 Hz)
# # QP_file = np.arange(0, 100, 1)     # 260 Hz
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
#
# QPT_Q6_Neg100_HMM = QPTunneling_Liu(name='Q6_Neg100_HMM')
# QPT_Q6_Neg100_HMM.add_datasets(QP_filenames, HMM=True)
#
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)



QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-05-20'
# experiment_name = ('Interleave_PSD_Neg16')
# # QP_file = np.arange(0, 10, 1)   # 233 H
# QP_file = np.arange(0, 100, 1)   # 232 Hz
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg16 = QPTunneling_Liu(name='Q6_Neg16')
# # QPT_Q6_Neg16 = QPTunneling(name='Q6_Neg16')
# QPT_Q6_Neg16.add_datasets(QP_filenames)

# date = '11-05-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 10, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames)
#
# date = '11-06-20'
# experiment_name = ('Interleave_PSD_Neg16')
# QP_file = np.arange(0, 10, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg16 = QPTunneling_Liu(name='Q6_Neg16')
# QPT_Q6_Neg16.add_datasets(QP_filenames)

# date = '11-19-20'
# experiment_name = ('Interleave_PSD_Neg100')
# # QP_file = np.arange(0, 200, 1)   #
# QP_file = np.arange(290, 300, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)

# date = '11-20-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(200, 400, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# date = '11-21-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# date = '11-21-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)

# date = '11-21-20'

# experiment_name = ('Old_Interleave_PSD_Neg100')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)
#
# experiment_name = ('Old_Interleave_PSD_Neg10')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# experiment_name = ('100us_Interleave_PSD_Neg100')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)
#
# experiment_name = ('100us_Interleave_PSD_Neg10')
# QP_file = np.arange(0, 130, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# date = '11-22-20'
#
# experiment_name = ('200us_Interleave_PSD_Neg100')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)
#
# experiment_name = ('200us_Interleave_PSD_Neg10')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# experiment_name = ('50us_Every10Cal_Interleave_PSD_Neg100')
# QP_file = np.arange(0, 20, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, simulate=True)
#
# experiment_name = ('50us_Every10Cal_Interleave_PSD_Neg10')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)

# date = '11-23-20'
#
# experiment_name = ('Interleave_PSD_Neg100_HiRO')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames, HMM=False)
# #
# experiment_name = ('Interleave_PSD_Neg10_HiRO')
# QP_file = np.arange(0, 200, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames, HMM=False)


date = '12-07-20'

# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg100 = QPTunneling_Wilen(name='Q6_Neg100')
# # QPT_Q6_Neg100 = QPTunneling_Liu(name='Q6_Neg100')
# QPT_Q6_Neg100.add_datasets(QP_filenames)
#
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg10 = QPTunneling_Wilen(name='Q6_Neg10')
# # QPT_Q6_Neg10 = QPTunneling_Liu(name='Q6_Neg10')
# QPT_Q6_Neg10.add_datasets(QP_filenames)
#
# date = '12-08-20'
# experiment_name = ('Interleave_PSD_Neg20')
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg20 = QPTunneling_Wilen(name='Q6_Neg20')
# QPT_Q6_Neg20.add_datasets(QP_filenames)
#
# date = '12-09-20'
# experiment_name = ('Interleave_PSD_Neg16')
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg16 = QPTunneling_Wilen(name='Q6_Neg16')
# QPT_Q6_Neg16.add_datasets(QP_filenames)
#
# experiment_name = ('Interleave_PSD_Neg13')
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_Neg13 = QPTunneling_Wilen(name='Q6_Neg13')
# QPT_Q6_Neg13.add_datasets(QP_filenames)



"""Some Q2 data"""

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q2_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-06-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 100, 1)   # 412 Hz
# # QP_file = np.arange(0, 10, 1)   # 412 Hz
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q2_Neg100 = QPTunneling_Liu(name='Q2_Neg100')
# QPT_Q2_Neg100.add_datasets(QP_filenames, simulate=False)

# date = '11-07-20'
# experiment_name = ('Interleave_PSD_Neg100')
# # QP_file = np.arange(0, 10, 1)   # 373 Hz
# # QP_file = np.arange(0, 50, 1)   # 321 Hz
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# # QPT_Q2_Neg100 = QPTunneling_Liu(name='Q2_Neg100')
# QPT_Q2_Neg100 = QPTunneling(name='Q2_Neg100')
# QPT_Q2_Neg100.add_datasets(QP_filenames)

# date = '11-07-20'
# experiment_name = ('Interleave_PSD_Neg20')
# # QP_file = np.arange(0, 10, 1)   # 343 Hz
# # QP_file = np.arange(0, 50, 1)   # 371 Hz
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q2_Neg20 = QPTunneling_Liu(name='Q2_Neg20')
# QPT_Q2_Neg20.add_datasets(QP_filenames)

# date = '11-08-20'
# experiment_name = ('Interleave_PSD_Neg16')
# # QP_file = np.arange(0, 50, 1)   # 391 Hz
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# # QPT_Q2_Neg16 = QPTunneling_Liu(name='Q2_Neg16')
# QPT_Q2_Neg16 = QPTunneling(name='Q2_Neg16')
# QPT_Q2_Neg16.add_datasets(QP_filenames)

# date = '11-08-20'
# experiment_name = ('Interleave_PSD_Neg13')
# # QP_file = np.arange(0, 50, 1)   # 340 Hz
# # QP_file = np.arange(50, 100, 1)   # 368 Hz
# QP_file = np.arange(0, 100, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# # QPT_Q2_Neg13 = QPTunneling_Liu(name='Q2_Neg13')
# QPT_Q2_Neg13 = QPTunneling(name='Q2_Neg13')
# QPT_Q2_Neg13.add_datasets(QP_filenames)

# date = '11-08-20'
# experiment_name = ('Interleave_PSD_Neg10')
# # QP_file = np.arange(0, 50, 1)   #  362 Hz
# # QP_file = np.arange(50, 100, 1)   # 378 Hz
# # QP_file = np.arange(0, 100, 1)   #
# QP_file = np.arange(0, 10, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q2_Neg10 = QPTunneling_Liu(name='Q2_Neg10')
# # QPT_Q2_Neg10 = QPTunneling(name='Q2_Neg10')
# QPT_Q2_Neg10.add_datasets(QP_filenames, simulate=False)

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q1_withQ2Poison/{}/{}/MATLABData/{}')
# date = '11-11-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 100, 1)   #434.68 Hz with HMM
# QP_file = np.arange(0, 10, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q1_Neg100 = QPTunneling_Liu(name='Q1_NoPoison')
# QPT_Q1_Neg100.add_datasets(QP_filenames, simulate=False, HMM=True)
# QPT_Q1_Neg100_HMM = QPTunneling_Liu(name='Q1_NoPoison_HMM')
# QPT_Q1_Neg100_HMM.add_datasets(QP_filenames, simulate=False, HMM=True)

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q1_withQ2Poison/{}/{}/MATLABData/{}')
# date = '11-15-20'
# experiment_name = ('Interleave_PSD_Neg20')
# QP_file = np.arange(0, 100, 1)   #
# # QP_file = np.arange(0, 10, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q1_Neg20 = QPTunneling_Liu(name='Q1_Neg20')
# QPT_Q1_Neg20.add_datasets(QP_filenames, HMM=True)

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q1_withQ2Poison/{}/{}/MATLABData/{}')
# date = '11-15-20'
# experiment_name = ('Interleave_PSD_Neg30')
# # QP_file = np.arange(0, 100, 1)   #
# # QP_file = np.arange(0, 10, 1)   #
# QP_file = np.arange(50, 60, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q1_Neg30 = QPTunneling_Liu(name='Q1_Neg30')
# QPT_Q1_Neg30.add_datasets(QP_filenames, HMM=True)

"""Q3 data"""
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q3_withQ2Poison/{}/{}/MATLABData/{}')
# experiment_name = ('50us_Interleave_PSD_Neg100')
# date = '12-15-20'
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# # QPT_Q6_Neg100 = QPTunneling_Wilen(name='Q6_Neg100')
# QPT_Q3_Neg100 = QPTunneling_Liu(name='Q3_Neg100_50us')
# QPT_Q3_Neg100.add_datasets(QP_filenames)

# experiment_name = ('25us_Interleave_PSD_Neg100')
# date = '12-16-20'
# QP_file = np.arange(0, 500, 1)   #
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q3_Neg100 = QPTunneling_Liu(name='Q3_Neg100_25us')
# QPT_Q3_Neg100.add_datasets(QP_filenames)

QPT_List = [QPT_Q4_Poison_Neg100dBm]
# QPT_List = [QPT_Q4_Poison_Neg100dBm_Clean, QPT_Q4_Poison_Neg100dBm_Dirty]
# QPT_List = [QPT_Q4_NoPoison, QPT_Q4_Poison, QPT_Q6_NoPoison]
# QPT_List = [QPT_Q6_NoPoison_Clean, QPT_Q6_NoPoison_Dirty]
# QPT_List = [QPT_Q2_Neg100, QPT_Q2_Neg20, QPT_Q2_Neg16, QPT_Q2_Neg13, QPT_Q2_Neg10]
# QPT_List = [QPT_Q4_NoPoison, QPT_Q4_Poison_Neg20dBm, QPT_Q4_Poison_Neg16dBm, QPT_Q4_Poison_Neg13dBm, QPT_Q4_Poison_Neg10dBm]
# QPT_List = [QPT_Q6_NoPoison, QPT_Q6_Neg20,  QPT_Q6_Neg16,  QPT_Q6_Neg13,  QPT_Q6_Neg10]
# QPT_List = [QPT_Q6_Neg100, QPT_Q6_Neg20, QPT_Q6_Neg16, QPT_Q6_Neg13, QPT_Q6_Neg10]
# QPT_List = [QPT_Q4_Poison]
# QPT_List = [QPT_Q2_Neg100]
# QPT_List = [QPT_Q2_Neg100, QPT_Q2_Neg20]
# QPT_List = [QPT_Q1_Neg100, QPT_Q1_Neg100_HMM]
# QPT_List = [QPT_Q4_Poison_Neg10dBm]
# QPT_List = [QPT_Q1_Neg30]
# QPT_List = [QPT_Q1_Neg100]
# QPT_List = [QPT_Q3_Neg100]

plotMultiFittedPSD(QPT_List)

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
# date = '10-30-20'
# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 200, 1)     # 363
# gamma_list = []
# for f in QP_file:
#     QPT = QPTunneling_Liu()
#     QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in np.arange(f, f+5, 1)]
#     QPT.add_datasets(QP_filenames)
#     QPT.get_psd(window_averaging=True)
#     QPT.get_fit_Liu()
#     gamma = QPT.params[0]
#     gamma_list.append(gamma)
#
# print(gamma_list)
#
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# ax.set_title('Gamma_vs_FileNumber_Neg10')
# ax.set_xlabel('File_Number')
# ax.set_ylabel('$S_\eta (\eta^2/Hz)$')
# ax.plot(gamma_list)
# plt.show()
