"""
2021Apr19
For 5 files P1 and 50 files of PSD
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty,QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np
import noiselib
import measurement.general.calibration as cal

#See if this works -- curve fitting with parameter uncertainties
import scipy.odr
import scipy.stats
from scipy.optimize import curve_fit


def getOccupation(IQData, IQCenter, IQCenterStd):
    WO = cal.WeightedOccupation()
    WO.importIQData(IQData)
    WO.importCalIQCenter(IQCenter)
    WO.importCalIQCenterStd(IQCenterStd)
    WO.plotIQData()

    WO.IQDataTo2DHistogram(bins=100)
    WO.fitIQDataHist()
    WO = WO.getOccupation()
    # print('WO=', WO)
    return WO


"""Q4"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/{}/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# First file -
user = 'LIU'
date = 'P1Pre2021May13'
temp = '425mK'

# Second file -
# user = 'LIU'
# date = '07-07-21'
# temp = '76mK'

experiment_name_P1 = ('Q4_P1_'+temp)
P1_file_Number = np.arange(0, 5, 1)
P1_file = [QP_path.format(user, date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in P1_file_Number]

# This gets me all the parameters--most importantly the IQ centers and stdevsp
variables = noiselib.loadmat_ExptVars(P1_file[0])
print(variables.keys())
# This gets me all the data
data=noiselib.loadmat(P1_file[0])
print(data.keys())

Is_pre=data['Is_pre']
Qs_pre=data['Qs_pre']

Ig_pre = variables['Ground_State_I_pre']
Qg_pre = variables['Ground_State_Q_pre']
Ie_pre = variables['Excited_State_I_pre']
Qe_pre = variables['Excited_State_Q_pre']

Ig_prestd = variables['Ground_State_I_prestd']
Qg_prestd = variables['Ground_State_Q_prestd']
Ie_prestd = variables['Excited_State_I_prestd']
Qe_prestd = variables['Excited_State_Q_prestd']

#this is not quite right--there are 36 sets of 5000 IQ points, each with a different Pre-Idle-RO time.
IQData_pre = np.column_stack((Is_pre, Qs_pre))
IQ_preCenter = np.array([[Ig_pre, Qg_pre], [Ie_pre, Qe_pre]])
IQ_preCenterStd = np.array([[Ig_prestd, Qg_prestd], [Ie_prestd, Qe_prestd]])

WeightedOccupation_pre = getOccupation(IQData=IQData_pre, IQCenter=IQ_preCenter, IQCenterStd=IQ_preCenterStd)


print(WeightedOccupation_pre)
print(len(IQData_pre))
