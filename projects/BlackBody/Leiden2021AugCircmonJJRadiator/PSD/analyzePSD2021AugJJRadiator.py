"""
2021Jul
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison, OneStateCleanDirty
import matplotlib.pyplot as plt
import numpy as np

"""Q2"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')
# Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\08-19-21
# date = '08-19-21'
# date = 'JJRadiatorQPT_2021Aug19'
date = 'JJRadiatorQPT_Q2Mode_2021Aug19'
QB_id = 'Q2'
# J2Biaslist = [40, 42, 44, 46, 48, 50, 55, 57, 59, 61, 63, 65, 67, 69, 71,
#               73, 75, 77, 79, 81, 83, 85]
# J2Q2Biaslist = np.arange(207, 250, 2)
J2Q2Biaslist = [85]
# J2Biaslist = [40, 42, 44]
J2_QPT_2D = []
# for J2Bias in J2Biaslist:
for J2Bias in J2Q2Biaslist:
    experiment_name_PSD = (QB_id+'_PSD_'+str(J2Bias)+'mDACJ2')
    PSD_file_Number = np.arange(0, 25, 1)
    PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
    QPT_Q = QPTunneling_Wilen(name='{} with J2 = {} mDAC, {}GHz'.
                              format(QB_id, str(J2Bias), str(J2Bias*4.8)))
    QPT_Q.add_datasets(PSD_file)

    QPT_List = [QPT_Q]
    plotMultiFittedPSD(QPT_List, save=False, name='{} with J2 ={} mDAC {}GHz'.
                              format(QB_id, str(J2Bias), str(J2Bias*4.8)))

    QPT_Q.get_fit()
    # print(QB_id)
    # print (J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
    J2_QPT_2D.append([J2Bias, QPT_Q.params[0]])

print('QB_id=', QB_id)
print('J2_QPT_2D=', J2_QPT_2D)