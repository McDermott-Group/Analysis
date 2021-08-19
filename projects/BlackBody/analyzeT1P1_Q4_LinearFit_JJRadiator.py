from antennalib import T1, P1, GammaUp
import numpy as np
import matplotlib.pyplot as plt

# Q3_data_array = Up_array()

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# date = 'JJRadiatorT1P1_2021Aug18'
date = 'JJRadiatorT1P1_2021Aug17'
experiment_name_P1_global = 'Q4_P1_'
experiment_name_T1_global = 'Q4_T1_'

GammaUp = GammaUp()

GammaUp_list=[]

data_type = 'Weighted_Occupation'

P1_avg = []
T1_avg = []

J2Bs = (list(np.arange(0,101,10)) + list(np.arange(59, 81, 2))) #August18
# J2Bs = list(np.arange(100,980,20)) #August17
J2Bs.sort()

for it,J2B in enumerate(J2Bs):
    print(it)
    experiment_name_P1 = experiment_name_P1_global+ str(J2B) + 'mDACJ2'
    experiment_name_T1 = experiment_name_T1_global + str(J2B) + 'mDACJ2'
    file_Number = np.arange(0, 5, 1)
    P1_file = [file_path.format(date, experiment_name_P1, experiment_name_P1) + '_{:03d}.mat'.format(i) for i in file_Number]
    P1_data = P1()
    P1_data.add_data_from_matlab(P1_file, 20,data_type1=data_type)
    P1_avg.append(np.average(P1_data.occ_1D_avg))
    T1_file = [file_path.format(date, experiment_name_T1, experiment_name_T1) + '_{:03d}.mat'.format(i) for i in file_Number]
    T1_data = T1()
    T1_data.add_data_from_matlab(T1_file, 20, data_type1=data_type, fit_type='Linear', num_points=5)
    # GammaUp.add_data(20, P1_data.occ_1D_avg, T1_data.Gamma_fit_parameters)
    # T1_data.plot()
    T1_avg.append(np.average(T1_data.Gamma_fit_parameters))
    # Gamma_Down_1D = GammaUp.Gamma_Down_2D[it]
    # P1_1D = GammaUp.P1_2D[it]
    # Gamma_Up_1D = P1_1D * Gamma_Down_1D / (1 - P1_1D)
    # Gamma_down_avg, Gamma_down_std = np.mean(Gamma_Down_1D), np.std(Gamma_Down_1D)
    # Gamma_up_avg, Gamma_up_std = np.mean(Gamma_Up_1D), np.std(Gamma_Up_1D)
    # P1_avg, P1_std = np.mean(P1_1D), np.std(P1_1D)
    # GammaUp_list.append(Gamma_up_avg)

plt.plot(J2Bs, P1_avg, 'o-', label='P1')
plt.xlabel('JJ Bias (mDAC)')
plt.ylabel('P1')
plt.grid()
plt.legend()
plt.show()

plt.plot(J2Bs, np.divide(1e6,T1_avg), 'o-', label='T1')
plt.xlabel('JJ Bias (mDAC)')
plt.ylabel('T1 (us)')
plt.grid()
plt.legend()
plt.show()

print(T1_avg)