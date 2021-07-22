"""
2021Jul
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q2"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'PSDP1T1Study2021Apr17'
date = '07-21-21'
QB_id = 'Q2'
# Poison = 'Poison_3dBm_20us'
Poison = 'Poison_Neg9dBm_20us'
# date = 'PSD2021Jun28'
# Q2_PSD_TestOne_000
# experiment_name_PSD = (QB_id+'_PSD_TestOne')
# experiment_name_PSD = (QB_id+'_PSD_TestTwo')
experiment_name_PSD = (QB_id+'_PSD_'+Poison)
PSD_file_Number = np.arange(0, 25, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
# QPT_Q4 = QPTunneling_Liu(name='{} at temp={} mK'.format(QB_id, temp))
QPT_Q2 = QPTunneling_Wilen(name='{} at {} Poison'.format(QB_id, Poison))
# QPT_Q2 = QPTunneling_Liu(name='{} at {} Poison'.format(QB_id, Poison))
QPT_Q2.add_datasets(PSD_file)


QPT_List = [QPT_Q2]
plotMultiFittedPSD(QPT_List)

QPT_Q2.get_fit()
print(QB_id)
print('temp (mK), QPT (Hz)')
# # # print (temp, QPT_Q4.params[0])
print ('[{}, {:.1f}]'.format(temp/1000.0, QPT_Q4.params[0]))