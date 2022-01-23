from cal_debuglib import blobCenterStd, TwoDFit_debug
import numpy as np
import matplotlib.pyplot as plt


file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'P1T12021Jun19'
date = '06-22-21'
experiment_name_P1_global = 'Q1_P1_'


# temp_list = [72, 79, 94, 103, 116, 124, 138, 145, 150, 160, 170, 178, 186,
#              193, 200, 205, 210, 216, 220, 225, 229, 234, 238]


temp_list = [80]
for temp in temp_list:
    Pe = []
    experiment_name_P1 = experiment_name_P1_global + str(temp) + 'mk'
    # for N in [16]:
    for N in range(25):
        file_Number = np.arange(N, N+1, 1)
        cal = TwoDFit_debug()
        cal_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
        cal.add_data_from_matlab(cal_file, temp)
        cal.fit1D_Hist()
        Pe.append(cal.Pe)
    print('Pe=', Pe)

iteration = np.arange(0, 25, 1)
plt.plot(iteration, Pe/np.mean(Pe), label='Pe')
plt.xlabel('Time')
plt.ylabel('Normalized')
plt.grid()
plt.legend()
plt.ylim([0.5, 2])
plt.show()

# plt.hist(Pe, 5)
# plt.show()












