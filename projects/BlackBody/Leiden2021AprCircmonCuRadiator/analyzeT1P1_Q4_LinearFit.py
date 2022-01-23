from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'P1T12021Aug05'
experiment_name_T1_global = 'Q4_T1_'
experiment_name_P1_global = 'Q4_P1_'
# experiment_name_T1_global = 'T1_'
# experiment_name_P1_global = 'P1_'
# experiment_name_T1_global = 'T1_Q1_Stats'


GammaUp = GammaUp()

temp_list_A = [80,103,126,144,151,167,182,200,214,230,252,273,299,332,365,398,445,500]
# temp_list=[238]
GammaUp_list=[]

data_type = 'Projected_Occupation'

for it,temp in enumerate(temp_list_A):
    print(it)
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mK'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mK'
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp, data_type1=data_type, fit_type='Linear', num_points=12)
    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in
               file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp, data_type1=data_type)
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

temp_list=temp_list_A

plt.plot(np.divide(1000.0,temp_list), np.log(np.divide(1.0,GammaUp_list)),'.', label='Gamma Up')
plt.xlabel('1/T (1/K)')
plt.ylabel(' Log(1/Gamma)')
plt.grid()
plt.legend()
plt.show()
