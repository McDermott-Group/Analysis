"""
2021Oct
Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\PSD2021Oct25_Q2WithJ7
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

QP_path = ('Z:/mcdermott-group/data/Antenna/SUXmon2/LIU/VitoChip2/{}/{}/MATLABData/{}')
# date = '10-20-21'
date = 'PSD2021Oct25_Q2WithJ7'
QB_id = 'Q1'
a1 = np.arange(0.0, 0.025, 0.005)
a2 = np.arange(0.025, 0.031, 0.0025)
# a3 = np.arange(0.04, 0.105, 0.005)
a = np.concatenate((a1, a2))
# JBlist = np.arange(34000, 64000, 2000)
# JBlist = [55000]
J_QPT_2D = []

for JB in a:
    experiment_name_PSD = (QB_id+'_PSD_'+str(int(1000000*JB))+'uDACJB')
    PSD_file_Number = np.arange(0, 25, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J7 = {} uDAC, {}GHz'.
                              format(QB_id, str(JB), str(JB*4.6)))
    # QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    QPT_Q.add_datasets(PSD_file)
    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=True, name='{} with J2 ={} uDAC {}GHz'.
                              format(QB_id, str(int(1000000*JB)), str(JB*4.6)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print(J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    rate = QPT_Q.params[0]  # units is Hz
    rate = int(rate*100)/100.0
    JB = int(JB*10000)/10000.0
    # print('rate=', rate)
    J_QPT_2D.append([JB, rate])
    # print('J_QPT_2D=', J_QPT_2D)

print('QB_id=', QB_id)
print('J_QPT_2D=', J_QPT_2D)


