# from projects.BlackBody.antennalib import T1
from antennalib import T1
import numpy as np

# Q3_data_array = Up_array()

QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = '06-18-21'
experiment_name_global = 'T1_Q1_Stats'

temp_list = [71]
for temp in temp_list:
    # experiment_name = experiment_name_global + str(temp) + 'mk'
    experiment_name = experiment_name_global

    file_Number_Tot = np.arange(0, 1, 1)
    file_Number_Extract = np.array([])
    file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
    #
    QP_file = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
    QP = T1()
    QP.add_data_from_matlab(QP_file, temp)
    # QP.plot_Chi()
    # QP.plot()
    # tempUpUpError = np.array([temp, QP.GammaUp, QP.GammaUp_std_err])
    # Q3_data_array.insert_data(tempUpUpError)
