from antennalib import QP_Up, Up_array
import numpy as np
import matplotlib.pyplot as plt

file_path = ('Z:/mcdermott-group/data/Antenna/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'J6Radiator_2021Oct01_UpRate'

# experiment_name = ('Q4_P1_54300uDACJ2')

def get_bias_list(qubit_number):
    a1 = np.arange(0.000, 0.016, 0.01)
    a2 = np.arange(0.016, 0.1, 0.004)
    a3 = np.arange(0.1, 0.15, 0.01)
    a = np.concatenate((a1, a2, a3))
    return a

data_dict = {}
for qubit_number in (1,2,4):
    data_dict[qubit_number] = []
    for J2B in get_bias_list(qubit_number):
        experiment_name =  'Q{}_P1_{}uDACJ2_T'.format(str(qubit_number), str(int(1000000 * J2B)))

        file_Number_Tt = np.arange(0, 5, 1)
        file_Number_Extract = np.array([])
        file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
        file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

        QP = QP_Up()
        QP.add_data_from_matlab(file_list, data_type1='Occupation_Filtered')
        #QP.plot_Chi()
        QP.plot(save=True,name='Figures/'+experiment_name+'.png')
        data_dict[qubit_number].append(QP.GammaUp)
    plt.plot(get_bias_list(qubit_number),data_dict[qubit_number],label='Q'+str(int(qubit_number)))
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('Gamma Up (Hz)')
plt.legend()
plt.show()

for qubit_number in (1,2,4):
    print get_bias_list(qubit_number)
    print data_dict[qubit_number]