from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'P1T12021Jun19'
date = '06-22-21'
experiment_name_T1_global = 'Q1_T1_'
experiment_name_P1_global = 'Q1_P1_'
# experiment_name_T1_global = 'T1_Q1_Stats'


GammaUp = GammaUp()
# temp_list = [72, 79, 94, 103, 116, 124, 138, 145, 150, 160, 170, 178, 186,
#              193, 200, 205, 210, 216, 220, 225, 229, 234, 238]

temp_list = [80]
for temp in temp_list:
    experiment_name_T1 = experiment_name_T1_global + str(temp) + 'mk'
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mk'
    file_Number = np.arange(0, 25, 1)
    #
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, temp)

    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, temp)
    # print('Gamma', T1_data.Gamma_fit_parameters)
    # print('P1', P1_data.occ_1D_avg)
    GammaUp.add_data(temp, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters)
    # QP.plot_Chi()
    # QP.plot()
    # QP.plot_histogram()
    # tempUpUpError = np.array([temp, QP.GammaUp, QP.GammaUp_std_err])

    # Q3_data_array.insert_data(tempUpUpError)

# print('temp', GammaUp.temp)
# print('P1_2D', GammaUp.P1_2D)
# print('Gamma_2D', GammaUp.Gamma_Down_2D)
# print('Gamma_Up_2D', GammaUp.Gamma_Up_2D)
# print('Gamma_Up_1D_error', GammaUp.Gamma_Up_1D_Error)

Gamma_Down_1D = GammaUp.Gamma_Down_2D[0]
P1_1D = GammaUp.P1_2D[0]
Gamma_Up_1D = P1_1D*Gamma_Down_1D/(1-P1_1D)
# print('Gamma_Down_1D=', Gamma_Down_1D)
Gamma_down_avg, Gamma_down_std = np.mean(Gamma_Down_1D), np.std(Gamma_Down_1D)
Gamma_up_avg, Gamma_up_std = np.mean(Gamma_Up_1D), np.std(Gamma_Up_1D)
P1_avg, P1_std = np.mean(P1_1D), np.std(P1_1D)
# print('Gamma_Up_1D=', Gamma_Up_1D)
print('Gamma_Down_1D=', Gamma_Down_1D)
print('P1_1D=', P1_1D)


# print('Gamma_up_avg=', Gamma_up_avg, 'Gamma_up_std=', Gamma_up_std)
# print('Gamma_down_avg=', Gamma_down_avg, 'Gamma_down_std=', Gamma_down_std)
# print('P1_avg=', P1_avg, 'P1_std=', P1_std)

iteration = np.arange(0, 25, 1)
plt.plot(iteration, Gamma_Down_1D/Gamma_down_avg, label='Gamma Down')
plt.plot(iteration, Gamma_Up_1D/Gamma_up_avg, label='Gamma Up')
plt.plot(iteration, P1_1D/P1_avg, label='P1')
plt.xlabel('Time')
plt.ylabel('Normalized')
plt.grid()
plt.legend()
plt.show()
