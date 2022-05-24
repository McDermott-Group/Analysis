from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'P1T12021Jun19'
date = '04-14-21'
# experiment_name_T1_global = 'Q1_T1_'
# experiment_name_P1_global = 'Q1_P1_'
experiment_name_T1_global = 'T1_'
experiment_name_P1_global = 'P1_'
# experiment_name_T1_global = 'T1_Q1_Stats'


GammaUp = GammaUp()

# temp_list = [72, 79, 94, 103, 116, 124, 138, 145, 150, 160, 170, 178, 186,
#              193, 200, 205, 210, 216, 220, 225, 229, 234, 238]

temp_list_A = [85,127,185,238,305,365,425,485,545]
# temp_list=[238]
GammaUp_list=[]

for it,temp in enumerate(temp_list_A):
    print(it)
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mK'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mK'
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp,fit_type='Linear',num_points=5)
    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    # print('Gamma', T1_data.Gamma_fit_parameters)
    # print('P1', P1_data.occ_1D_avg)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters)
    T1_data.plot()
    Gamma_Down_1D = GammaUp.Gamma_Down_2D[it]
    P1_1D = GammaUp.P1_2D[it]
    Gamma_Up_1D = P1_1D * Gamma_Down_1D / (1 - P1_1D)
    # print('Gamma_Down_1D=', Gamma_Down_1D)
    Gamma_down_avg, Gamma_down_std = np.mean(Gamma_Down_1D), np.std(Gamma_Down_1D)
    Gamma_up_avg, Gamma_up_std = np.mean(Gamma_Up_1D), np.std(Gamma_Up_1D)
    P1_avg, P1_std = np.mean(P1_1D), np.std(P1_1D)
    # print('Gamma_Up_1D=', Gamma_Up_1D)
    GammaUp_list.append(Gamma_up_avg)


temp_list_B = [76,82,95,105,110,116,125,139,145,155,163,170,175,180,186]
# temp_list=[238]
date = 'PSDP1T1Study2021Apr21'
# experiment_name_T1_global = 'Q1_T1_'
# experiment_name_P1_global = 'Q1_P1_'
experiment_name_T1_global = 'Interleave_T1_'
experiment_name_P1_global = 'Interleave_P1_'

offset= len(GammaUp_list)
for it,temp in enumerate(temp_list_B):
    print(it+offset)
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mK'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mK'
    file_Number = np.arange(0, 5, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp,fit_type='Linear',num_points=5)
    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    # print('Gamma', T1_data.Gamma_fit_parameters)
    # print('P1', P1_data.occ_1D_avg)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters)
    T1_data.plot()
    Gamma_Down_1D = GammaUp.Gamma_Down_2D[it+offset]
    P1_1D = GammaUp.P1_2D[it+offset]
    Gamma_Up_1D = P1_1D * Gamma_Down_1D / (1 - P1_1D)
    # print('Gamma_Down_1D=', Gamma_Down_1D)
    Gamma_down_avg, Gamma_down_std = np.mean(Gamma_Down_1D), np.std(Gamma_Down_1D)
    Gamma_up_avg, Gamma_up_std = np.mean(Gamma_Up_1D), np.std(Gamma_Up_1D)
    P1_avg, P1_std = np.mean(P1_1D), np.std(P1_1D)
    # print('Gamma_Up_1D=', Gamma_Up_1D)
    GammaUp_list.append(Gamma_up_avg)


temp_list_C = [80,185,200,210,220,229,237,245,252,260,266,274,284,294,305,318,331,342,353,364,375,386,396,404,412,424,438,449,462,472,484]
# temp_list=[238]
date = 'PSDP1T1Study2021Apr17'
# experiment_name_T1_global = 'Q1_T1_'
# experiment_name_P1_global = 'Q1_P1_'
experiment_name_T1_global = 'Interleave_T1_'
experiment_name_P1_global = 'Interleave_P1_'

offset= len(GammaUp_list)
for it,temp in enumerate(temp_list_C):
    print(it+offset)
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mK'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mK'
    file_Number = np.arange(0, 5, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp,fit_type='Linear',num_points=5)
    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    # print('Gamma', T1_data.Gamma_fit_parameters)
    # print('P1', P1_data.occ_1D_avg)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters)
    T1_data.plot()
    Gamma_Down_1D = GammaUp.Gamma_Down_2D[it+offset]
    P1_1D = GammaUp.P1_2D[it+offset]
    Gamma_Up_1D = P1_1D * Gamma_Down_1D / (1 - P1_1D)
    # print('Gamma_Down_1D=', Gamma_Down_1D)
    Gamma_down_avg, Gamma_down_std = np.mean(Gamma_Down_1D), np.std(Gamma_Down_1D)
    Gamma_up_avg, Gamma_up_std = np.mean(Gamma_Up_1D), np.std(Gamma_Up_1D)
    P1_avg, P1_std = np.mean(P1_1D), np.std(P1_1D)
    # print('Gamma_Up_1D=', Gamma_Up_1D)
    GammaUp_list.append(Gamma_up_avg)

temp_list=temp_list_A+temp_list_B+temp_list_C

plt.plot(np.divide(1000.0,temp_list), np.log(np.divide(1.0,GammaUp_list)),'.', label='Gamma Up')
plt.xlabel('1/T (1/K)')
plt.ylabel(' Log(1/Gamma)')
plt.grid()
plt.legend()
plt.show()

for temp in temp_list:
    print(temp)

for rate in GammaUp_list:
    print(np.log(np.divide(1.0,rate)))