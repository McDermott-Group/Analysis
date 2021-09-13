"""
2021June
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

"""Q4"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/DaveHarrison/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'PSDP1T1Study2021Apr17'
date = '07-12-21'
QB_id = 'Q2'
# date = 'PSD2021Jun28'
temp = 76

experiment_name_PSD = (QB_id+'_PSD_'+str(temp)+'mK_retunedQubitFreq')
# PSD_file_Number = np.arange(89, 189, 1)
PSD_file_Number = np.arange(0, 50, 1)
PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
QPT_Q4 = QPTunneling_Harrison(name='{} at temp={} mK'.format(QB_id, temp))
# QPT_Q4 = QPTunneling_Wilen(name='{} at temp={} mK'.format(QB_id, temp))
QPT_Q4.add_datasets(PSD_file)



QPT_List = [QPT_Q4]
# QPT_List = [QPT_Q4_temp2, QPT_Q4_temp3, QPT_Q4_temp4, QPT_Q4_temp5]
plotMultiFittedPSD(QPT_List)

QPT_Q4.get_fit()
print(QB_id)
print('temp (mK), QPT (Hz)')
# # # print (temp, QPT_Q4.params[0])
print ('[{}, {:.1f}]'.format(temp/1000.0, QPT_Q4.params[0]))