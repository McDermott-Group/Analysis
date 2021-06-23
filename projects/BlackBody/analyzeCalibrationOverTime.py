from cal_debuglib import blobCenterStd, TwoDFit_debug
import numpy as np
import matplotlib.pyplot as plt


file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'P1T12021Jun19'
experiment_name_P1_global = 'Q1_P1_'


# temp_list = [72, 79, 94, 103, 116, 124, 138, 145, 150, 160, 170, 178, 186,
#              193, 200, 205, 210, 216, 220, 225, 229, 234, 238]

N = 25

temp_list = [178]
for temp in temp_list:
    cal = blobCenterStd()
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mk'
    file_Number = np.arange(0, N, 1)

    cal_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    cal.add_data_from_matlab(cal_file, temp)

# print(cal.eIavg)
# plt.figure(1)
# iteration = np.arange(0, N, 1)
# plt.plot(iteration, cal.gIavg, label='g_I')
# plt.plot(iteration, cal.gQavg, label='g_Q')
# plt.plot(iteration, cal.eIavg, label='e_I')
# plt.plot(iteration, cal.eQavg, label='e_Q')
# plt.fill_between(iteration, cal.gIavg-cal.gIstd*0.1, cal.gIavg+cal.gIstd*0.1, alpha=0.1)
# plt.fill_between(iteration, cal.gQavg-cal.gQstd*0.1, cal.gQavg+cal.gQstd*0.1, alpha=0.1)
# plt.fill_between(iteration, cal.eIavg-cal.eIstd*0.1, cal.eIavg+cal.eIstd*0.1, alpha=0.1)
# plt.fill_between(iteration, cal.eQavg-cal.eQstd*0.1, cal.eQavg+cal.eQstd*0.1, alpha=0.1)
# plt.xlabel('Time')
# plt.ylabel('V +- stand deviation*0.1')
# plt.grid()
# plt.legend()
# # plt.show()
#
# """IQ"""
# plt.figure(2)
# plt.scatter(cal.gIavg, cal.gQavg, label='ground')
# plt.scatter(cal.eIavg, cal.eQavg, label='excited')
# plt.xlabel('I')
# plt.ylabel('Q')
# plt.grid()
# plt.legend()
# plt.show()