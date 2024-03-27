from antennalib import QP_Up, Up_array
import numpy as np
import matplotlib.pyplot as plt

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorP1Pre_2021Sep01'

# experiment_name = ('Q4_P1_54300uDACJ2')

def get_bias_list(qubit_number):
    if qubit_number ==1:
        return np.array([
            0,10,16.3,18.3,20.3,22.3,24.3,26.3,28.3,
            30.3,32.299,34.3,36.3,38.3,40.3,42.3,44.3,
            46.3,48.3,50.3,52.3,54.3,56.3,58.3,60.3,
            75,95,115,135])
    elif qubit_number == 2:
        return np.array([0,10]+list(np.arange(16.3,60.301,2))+list(np.arange(69.8,229.801,8)))
    elif qubit_number == 4:
        return [32.3,28.3]
        #return np.array([0,10]+list(np.arange(16.3,60.301,2))+list(np.arange(75,135.001,20)))

data_dict = {}
error_dict = {}
for qubit_number in [4]:#(1,2,4):
    data_dict[qubit_number] = []
    error_dict[qubit_number] = []
    for J2B in get_bias_list(qubit_number):
        experiment_name = ('Q'+str(qubit_number)+'_P1_'+str(int(J2B*1000))+'uDACJ2')

        file_Number_Tot = np.arange(0, 20, 1)
        file_Number_Extract = np.array([])
        file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
        file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

        QP = QP_Up()
        QP.add_data_from_matlab(file_list, data_type1='Occupation_Filtered')
        #QP.plot_Chi()
        QP.plotErrorBars(save=True,name='Figures/NEW'+experiment_name+'.png')
        data_dict[qubit_number].append(QP.GammaUp)
        error_dict[qubit_number].append(QP.GammaUp_std_err)
    plt.plot(get_bias_list(qubit_number),data_dict[qubit_number],label='Q'+str(int(qubit_number)))
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('Gamma Up (Hz)')
plt.legend()
plt.show()

for qubit_number in (1,2,4):
    print(get_bias_list(qubit_number))
    print(data_dict[qubit_number])
    print(error_dict[qubit_number])