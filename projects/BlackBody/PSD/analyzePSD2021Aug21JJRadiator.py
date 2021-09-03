"""
2021Jul
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q2"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
date = 'JJRadiatorQPT_Q1_AroundPeakFast_2021Aug29'
QB_id = 'Q1'
J2Q2Biaslist = np.arange(20000, 38000, 500)

J2_QPT_2D = []
for J2Bias in J2Q2Biaslist:
    experiment_name_PSD = (QB_id+'_PSD_'+str(J2Bias)+'uDACJ2')
    PSD_file_Number = np.arange(0, 50, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J2 = {} uDAC, {}GHz'.
                              format(QB_id, str(J2Bias), str(J2Bias*4.8)))
    QPT_Q.add_datasets(PSD_file)

    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=True, name='{} with J2 ={} mDAC {}GHz'.
                              format(QB_id, str(J2Bias), str(J2Bias*4.8)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print (J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    J2_QPT_2D.append([J2Bias, QPT_Q.params[0]])

print('QB_id=', QB_id)
print('J2_QPT_2D=', J2_QPT_2D)