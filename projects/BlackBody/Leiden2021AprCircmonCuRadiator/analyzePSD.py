"""
2021June
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'PSDP1T1Study2021Apr17'
# date = '06-27-21'
QB_id = 'Q4'
date = 'PSD2021Jun28'
# temp = 3
# temp = 76
# temp = 102
# temp = 125 #Q3
# temp = 199
# temp = 251
temp = 284
# temp = 331
# temp = 372
# temp = 404
# temp = 451
# temp = 501

experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK')
if QB_id == 'Q1' and temp == '76':
    PSD_file_Number = np.arange(100, 200, 1)
elif QB_id == 'Q3' and temp == '76':
    PSD_file_Number = np.arange(200, 400, 1)
else:
    PSD_file_Number = np.arange(0, 100, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4 = QPTunneling_Liu(name='{} at temp={} mK'.format(QB_id, temp))
QPT_Q4 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
QPT_Q4.add_datasets(PSD_file)

# QPT_Q4_temp1 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4_temp1.add_datasets(PSD_file)

# temp = 76
# experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK')
# PSD_file_Number = np.arange(200, 400, 1)
# PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4_temp2 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4_temp2.add_datasets(PSD_file)
# #
# temp = 125
# experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK')
# PSD_file_Number = np.arange(0, 100, 1)
# PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4_temp3 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4_temp3.add_datasets(PSD_file)
#
# temp = 199
# experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK')
# PSD_file_Number = np.arange(0, 100, 1)
# PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4_temp4 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4_temp4.add_datasets(PSD_file)
#
# temp = 284
# experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK')
# PSD_file_Number = np.arange(0, 100, 1)
# PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4_temp5 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4_temp5.add_datasets(PSD_file)


QPT_List = [QPT_Q4]
# QPT_List = [QPT_Q4_temp2, QPT_Q4_temp3, QPT_Q4_temp4, QPT_Q4_temp5]
plotMultiFittedPSD(QPT_List)

QPT_Q4.get_fit()
print(QB_id)
print('temp (mK), QPT (Hz)')
# # # print(temp, QPT_Q4.params[0])
print('[{}, {:.1f}]'.format(temp/1000.0, QPT_Q4.params[0]))