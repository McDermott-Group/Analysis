from SFQlib import T1_QP_1D, T2_1D, T1_QP_2D
import numpy as np


# """
# T1 1D Fit
# """
# file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# date = '12-28-21'
# experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep_1D'
# # experiment_name_T1 = 'T1'
#
# # file_Number = np.arange(20, 21, 1)
# # file_Number = np.arange(2, 3, 1)
# file_Number = [1]
# T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
# T1_data = T1_QP_1D()
# # T1_data.add_data_from_matlab(T1_file, fit_type='QP')
# T1_data.add_data_from_matlab(T1_file)
# T1_data.plot()

"""
T1 2D Fit
"""
# file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# # date = '12-28-21'
# date = '05-17-22'
# experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep'
# # experiment_name_T1 = 'T1_SFQ_Poison_Recovery'
# file_Number = [2]
# T1_2D_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
# T1_2D_data = T1_QP_2D()
# T1_2D_data.add_data_from_matlab(T1_2D_file)
# T1_data.plot()

"""
T2 fit
"""
file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# date = '12-29-21'
date = '05-19-22'
# experiment_name_T2 = 'T2_SFQ_Poison_Time_Sweep_1D_Study'
# experiment_name_T2 = 'T2_SFQ_Poison_During_Study'
experiment_name_T2 = 'T2_SFQ_Poison_During_Study_Q2'
# file_Number = np.arange(20, 21, 1)
file_Number = np.arange(0, 5, 1)
# file_Number = np.arange(0, 1, 1)
# file_Number = [30]
T2_file = [file_path.format(date, experiment_name_T2, experiment_name_T2) + '_{:03d}.mat'.format(i) for i in file_Number]
T2_data = T2_1D()
T2_data.add_data_from_matlab(T2_file)
# T2_data.plot()
fitParameters = T2_data.fit_parameters
T2_data.plot_Detuning()
# T2_data.plot()
# print(fitParameters)


