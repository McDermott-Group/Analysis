from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorT1_2021Aug26_HighDensity'
experiment_name_Q1_T1_global = 'Q1_T1_'
experiment_name_Q2_T1_global = 'Q2_T1_'
experiment_name_Q3_T1_global = 'Q3_T1_'
experiment_name_Q4_T1_global = 'Q4_T1_'

file_Number = np.arange(0, 50, 1)

data_type = 'Projected_Occupation'

T1_J2_Q1 =[]
T1_J2_Q2 =[]
T1_J2_Q3 =[]
T1_J2_Q4 =[]

J2_Bias_Q1=sorted(list(np.arange(0,20,10)*1000)+list(np.arange(20,40,1)*1000)+list(np.arange(40,100,10)*1000))
suffix = 'uDACJ2'
for it,J2B in enumerate(J2_Bias_Q1):
    print(it)
    experiment_name_T1 = experiment_name_Q1_T1_global + str(J2B) + suffix
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 20, data_type1=data_type, fit_type='Linear', num_points=20)
    T1_J2_Q1.append(1e6/np.average(T1_data.Gamma_fit_parameters))

J2_Bias_Q2=list(np.arange(0,75,10)*1000)+list(np.arange(75,100,1)*1000)
J2_Bias_Q2.sort()
suffix = 'uDACJ2'
for it,J2B in enumerate(J2_Bias_Q2):
    print(it)
    experiment_name_T1 = experiment_name_Q2_T1_global + str(J2B) + suffix
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 20, data_type1=data_type, fit_type='Linear', num_points=20)
    T1_J2_Q2.append(1e6/np.average(T1_data.Gamma_fit_parameters))

J2_Bias_Q3=list(np.arange(0,6,5)*1000)+list(np.arange(5,25,1)*1000)+list(np.arange(25,100,10)*1000)
J2_Bias_Q3.sort()
suffix = 'uDACJ2'
for it,J2B in enumerate(J2_Bias_Q3):
    print(it)
    experiment_name_T1 = experiment_name_Q3_T1_global + str(J2B) + suffix
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 20, data_type1=data_type, fit_type='Linear', num_points=20)
    T1_J2_Q3.append(1e6/np.average(T1_data.Gamma_fit_parameters))

J2_Bias_Q4=list(np.arange(0,40,10)*1000)+list(np.arange(40,60,1)*1000)+list(np.arange(60,100,10)*1000)
J2_Bias_Q4.sort()
suffix = 'uDACJ2'
for it,J2B in enumerate(J2_Bias_Q4):
    print(it)
    experiment_name_T1 = experiment_name_Q4_T1_global + str(J2B) + suffix
    file_Number = np.arange(0, 10, 1)
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 20, data_type1=data_type, fit_type='Linear', num_points=20)
    T1_J2_Q4.append(1e6/np.average(T1_data.Gamma_fit_parameters))



print(len(J2_Bias_Q1))
print(len(T1_J2_Q1))
print(len(J2_Bias_Q2))
print(len(T1_J2_Q2))
print(len(J2_Bias_Q3))
print(len(T1_J2_Q3))
print(len(J2_Bias_Q4))
print(len(T1_J2_Q4))

J2_Bias_Q1[:] = [x / 1000 for x in J2_Bias_Q1]
J2_Bias_Q2[:] = [x / 1000 for x in J2_Bias_Q2]
J2_Bias_Q3[:] = [x / 1000 for x in J2_Bias_Q3]
J2_Bias_Q4[:] = [x / 1000 for x in J2_Bias_Q4]

print('J2_Bias_Q1=', J2_Bias_Q1)
print('J2_Bias_Q2=', J2_Bias_Q2)
print('J2_Bias_Q3=', J2_Bias_Q3)
print('J2_Bias_Q4=', J2_Bias_Q4)
print('T1_J2_Q1=', T1_J2_Q1)
print('T1_J2_Q2=', T1_J2_Q2)
print('T1_J2_Q3=', T1_J2_Q3)
print('T1_J2_Q4=', T1_J2_Q4)


plt.plot(J2_Bias_Q1, T1_J2_Q1, label='Q1')
plt.plot(J2_Bias_Q2, T1_J2_Q2, label='Q2')
# plt.plot(J2_Bias_Q3, T1_J2_Q3, label='Q3')
plt.plot(J2_Bias_Q4, T1_J2_Q4, label='Q4')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('T1 (us)')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()


