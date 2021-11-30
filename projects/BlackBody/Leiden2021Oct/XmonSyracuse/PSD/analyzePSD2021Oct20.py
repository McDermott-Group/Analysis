"""
2021Oct
Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\10-20-21
Z:\mcdermott-group\data\Antenna\SUXmon\Liu\VitoChip1\PSD2021Oct20_Q123WithJ7
Z:\mcdermott-group\data\Antenna_Temporary\SUXmon\Liu\VitoChip1\11-01-21\Q2_PSD_0uDACJB_Transient\MATLABData
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

QP_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')
# date = '10-29-21'
date = '11-20-21'
# date = 'PSD2021Oct20_Q123WithJ7'
QB_id = 'Q1'
# a1 = np.arange(0.0, 0.025, 0.005)
# a2 = np.arange(0.025, 0.04, 0.0025)
# a3 = np.arange(0.04, 0.07, 0.001)
# a4 = np.arange(0.07, 0.1, 0.0025)
# a5 = np.arange(0.1, 0.15, 0.0025)
# a6 = np.arange(0.15, 0.25, 0.005)
# a = np.concatenate((a1, a2, a3, a4, a5, a6))
# JBlist = np.arange(34000, 64000, 2000)
# JBlist = [55000]
a = [0.100]
J_QPT_2D = []

for JB in a:
    experiment_name_PSD = (QB_id+'_PSD_'+str(int(1000000*JB))+'uDACJ1')
    PSD_file_Number = np.arange(10, 15, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J1 = {} uDAC, {}GHz'.
                              format(QB_id, str(JB), str(JB*4.6)))
    # QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    QPT_Q.add_datasets(PSD_file)
    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=False, name='{} with J2 ={} uDAC {}GHz'.
                              format(QB_id, str(int(1000000*JB)), str(JB*4.6)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print (J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    rate = QPT_Q.params[0]  # units is Hz
    rate = int(rate*100)/100.0
    JB = int(JB*10000)/10000.0
    # print('rate=', rate)
    J_QPT_2D.append([JB, rate])
    # print('J_QPT_2D=', J_QPT_2D)

print('QB_id=', QB_id)
print('J_QPT_2D=', J_QPT_2D)


