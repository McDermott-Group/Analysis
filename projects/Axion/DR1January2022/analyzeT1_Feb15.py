from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/testProject/Keysight/DCH/NA/{}/{}/MATLABData/{}')

date = '02-15-22'
experiment_name_Q3_T1_global = 'T1_Q3_'

data_type = 'Projected_Occupation'

T1_J2_Q3 =[]

J2_Bias_Q3=list(np.arange(0,486,2))
J2_Bias_Q3.sort()
suffix = 'mVJ2'
for it,J2B in enumerate(J2_Bias_Q3):
    print(it)
    experiment_name_T1 = experiment_name_Q3_T1_global + str(J2B) + suffix
    print(experiment_name_T1)
    file_Number = np.arange(0, 1, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 40, data_type1=data_type, num_points=40)
    T1_J2_Q3.append(1e6/np.average(T1_data.Gamma_fit_parameters))
    print((T1_data.Gamma_fit_parameters))

J2_Bias_Q3[:] = [484*2*(x-30)/1000 for x in J2_Bias_Q3]

plt.plot(J2_Bias_Q3, T1_J2_Q3, label='Q3')
plt.xlabel('J2 Freq (GHz)')
plt.ylabel('T1 (us)')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()


