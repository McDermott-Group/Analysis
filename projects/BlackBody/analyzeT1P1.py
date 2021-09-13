from antennalib import T1, P1, GammaUp, GammaUp_New
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'P1T12021Jun23'
# date = '06-22-21'
QB_id = 'Q4'
experiment_name_T1_global = QB_id+'_T1_'
experiment_name_P1_global = QB_id+'_P1_'
# experiment_name_T1_global = 'T1_Q1_Stats'


GammaUp = GammaUp_New()
# GammaUp = GammaUp()

temp_list = [73, 124, 142, 167, 199, 215, 222, 230, 237, 242]
# temp_list = [222]
for temp in temp_list:
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mk'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mk'
    file_Number = np.arange(0, 40, 1)
    #
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp)

    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters, QB_id)
# GammaUp.plot_iteration_at_temp()
date = 'P1T12021Jun26'
temp_list = [92, 102, 111, 252, 332]
for temp in temp_list:
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mk'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mk'
    file_Number = np.arange(0, 40, 1)
    #
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp)

    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters, QB_id)

GammaUp.plot()