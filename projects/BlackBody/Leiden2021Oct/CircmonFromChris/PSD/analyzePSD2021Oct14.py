"""
2021Oct
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np
QP_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
# date = '10-15-21'
# date = '2021Oct16_QB124_PSD_J6Radiator'
date = '2021Oct19_QB2_PSD_J6Radiator_HiFreq_10usRep'
# date = '2021Oct19_QB4_PSD_RepRate'
QB_id = 'Q2'
# a1 = np.arange(83000, 100000, 1000)
# a2 = np.arange(100000, 150000, 2500)
# JBlist = np.concatenate((a1, a2))
# JBlist = np.arange(350000, 450000, 20000)
# JBlist = [150000, 169999, 189999, 209999, 229999, 249999]
# JBlist = [269999, 289999, 309999, 329999]
JBlist = [270000]
J_QPT_2D = []

for JB in JBlist:
    experiment_name_PSD = (QB_id+'_PSD_'+str(JB)+'uDACJ6')
    PSD_file_Number = np.arange(0, 50, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J6 = {} uDAC, {}GHz'.
                              format(QB_id, str(JB), str(JB*4.6)))
    # QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    QPT_Q.add_datasets(PSD_file)
    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=True, name='{} with J2 ={} uDAC {}GHz'.
                              format(QB_id, str(JB), str(JB*4.6)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print (J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    rate = QPT_Q.params[0]  # units is Hz
    rate = int(rate*100)/100.0
    # print('rate=', rate)
    J_QPT_2D.append([JB, rate])
    # print('J_QPT_2D=', J_QPT_2D)

print('QB_id=', QB_id)
print('J_QPT_2D=', J_QPT_2D)

