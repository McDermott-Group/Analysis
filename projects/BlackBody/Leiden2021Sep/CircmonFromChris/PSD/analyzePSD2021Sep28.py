"""
2021Sep
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np

QP_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
# date = 'J6Radiator_2021Sep28_Q124'
date = 'J6Radiator_2021Sep28_Q2_Retake_2'
QB_id = 'Q2'
# a1 = np.arange(0, 16000, 5000)
# a2 = np.arange(16000, 44000, 2000)
# JBlist = np.concatenate((a1, a2))
# JBlist = np.arange(34000, 64000, 2000)
JBlist = [240000]
J_QPT_2D = []

for JB in JBlist:
    experiment_name_PSD = (QB_id+'_PSD_'+str(JB)+'uDACJB')
    PSD_file_Number = np.arange(20, 25, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J6 = {} uDAC, {}GHz'.
                              format(QB_id, str(JB), str(JB*4.6)))
    # QPT_Q = QPTunneling_Liu(name='{} with J6 = {} uDAC, {}GHz'.
    #                           format(QB_id, str(JB), str(JB*4.6)))
    QPT_Q.add_datasets(PSD_file)
    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=False, name='{} with J2 ={} uDAC {}GHz'.
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



"""
[0, 371.14], [2000, 369.97], [4000, 375.69], [6000, 356.82], [8000, 369.57]
, [10000, 353.39], [12000, 365.89], [14000, 609.26], [16000, 1071.49], [18000, 2324.16], [
20000, 2266.18], [22000, 4337.28], [24000, 7424.02], [26000, 5392.25], [28000, 2161.45], [
30000, 1818.08], [32000, 1081.65], [34000, 684.81], [36000, 691.99], [38000, 759.46], [400
00, 898.52], [42000, 1081.35], [44000, 1354.85], [46000, 1735.18], [48000, 1948.75], [5000
0, 2154.12], [52000, 2497.86], [54000, 4054.68], [56000, 4720.82], [58000, 5113.44], [6000
0, 6098.86], [62000, 8058.75], [64000, 9301.08], [66000, 6947.42], [68000, 9125.52], [7000
0, 7989.44], [72000, 13593.55], [74000, 18918.19], [76000, 10437.61], [78000, 16761.9], [8
0000, 37059.05], [82000, 40321.95], [84000, 40401.47], [86000, 43726.85], [88000, 54073.7]
, [90000, 41859.17], [92000, 37337.43], [94000, 52313.14], [96000, 44321.85], [98000, 5234
7.23]]
"""