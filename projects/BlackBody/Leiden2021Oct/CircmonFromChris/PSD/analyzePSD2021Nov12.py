"""
2021Oct
Z:\mcdermott-group\data\Antenna\Circmon\Liu\CW20180514A_Ox2\2021Nov12_QB12_PSD_J6Radiator_Sim
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np
QP_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
# date = '11-08-21'
# date = '2021Nov08_Q3_PSD_J6Radiator'
# date = '2021Nov12_QB12_PSD_J6Radiator_Sim'
# date = '2021Nov12_QB14_PSD_Flux_J6Radiator_Sim'
# date = '2021Nov15_QB4_PSD_J6Radiator_Sim'
date = '2021Nov16_QB4_PSD_J6Radiator_Sim'
# date = '2021Nov02_QB124_PSD_J6Radiator_Sim'
# date = '2021Oct19_QB2_PSD_J6Radiator_HiFreq_10usRep'
# date = '2021Oct19_QB4_PSD_RepRate'
QB_id = 'Q2'
# a1 = np.arange(83000, 100000, 1000)
# a2 = np.arange(100000, 150000, 2500)
# JBlist = np.concatenate((a1, a2))
# JBlist = np.arange(350000, 450000, 20000)
# JBlist = [150000, 169999, 189999, 209999, 229999, 249999]
# JBlist = [269999, 289999, 309999, 329999]
# a1 = np.arange(0.0, 0.1, 0.010)
# a2 = np.arange(0.1, 0.5, 0.005)
# a = np.concatenate((a1, a2))
# a = np.arange(0.300, 0.299, 0.001)
# a1 = np.arange(0.300, 0.340, 0.001)
# a2 = np.arange(0.6, 0.65, 0.001)
# a = np.concatenate((a1, a2))
# a = np.arange(0.5, 0.7, 0.005)
# a = np.arange(0.511, 0.530, 0.001)
# a = np.arange(0.0, 0.5, 0.01)
# a = -np.arange(0.24, 0.285, 0.01)
# a1 = np.arange(0.100, 0.190, 0.001)
a = np.arange(0.190, 0.260, 0.001)
# a = np.arange(0.605, 0.7, 0.005)
# a = np.concatenate((a1, a2))
# a = a2
JBlist = [int(j*1000) for j in a]
# JBlist = [0]
J_QPT_2D = []

for JB in JBlist:
    experiment_name_PSD = (QB_id+'_PSD_'+str(JB)+'mVJ6_HighDensity')
    # experiment_name_PSD = (QB_id+'_PSD_Q1Q4LSS')
    PSD_file_Number = np.arange(0, 25, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J6 = {} mV, {}GHz'.
                              format(QB_id, str(JB), str(JB*0.97)))
    # QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    # QPT_Q = QPTunneling_Wilen(name='{}_USS'.format(QB_id))
    QPT_Q.add_datasets(PSD_file)
    QPT_List = [QPT_Q]
    # plotMultiFittedPSD(QPT_List, save=True, name='{} with J2 ={} uDAC {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    plotMultiFittedPSD(QPT_List, save=True, name='{} with J6 = {} mV, {}GHz'.
                              format(QB_id, str(JB), str(JB*0.97)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print(J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    rate = QPT_Q.params[0]  # units is Hz
    rate = int(rate*100)/100.0
    # print('rate=', rate)
    J_QPT_2D.append([JB, rate])
    # print('J_QPT_2D=', J_QPT_2D)

print('QB_id=', QB_id)
print('J_QPT_2D=', J_QPT_2D)


