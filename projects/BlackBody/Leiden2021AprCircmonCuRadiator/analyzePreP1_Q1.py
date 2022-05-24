from antennalib import QP_Up, Up_array
import numpy as np

Q3_data_array = Up_array()

QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')


date = 'P1Pre2021Jun11'
# date = 'P1Pre2021Jun18'
experiment_name_global = 'Q3_P1_'

# temp_list = [78, 118, 126, 138, 151, 170, 194, 200, 237, 274, 306, 336, 365, 397, 426]
temp_list = [126]
for temp in temp_list:
    experiment_name = experiment_name_global + str(temp) + 'mk'

    file_Number_Tot = np.arange(0, 50, 1)
    file_Number_Extract = np.array([])
    file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
    #
    QP_file = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
    QP = QP_Up()
    QP.add_data_from_matlab(QP_file, temp)
    QP.plot_Chi()
    QP.plot()
    # tempUpUpError = np.array([temp, QP.GammaUp, QP.GammaUp_std_err])
    # Q3_data_array.insert_data(tempUpUpError)


# experiment_name = experiment_name_global + '185mk_new'
# temp = '185mk_new'
# file_Number_Tot = np.arange(0, 5, 1)
# file_Number_Extract = np.array([])
# file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
# #
# QP_file = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
# QP = QP_Up()
# QP.add_data_from_matlab(QP_file, temp)
# # QP.plot()
# tempUpUpError = np.array([185, QP.GammaUp, QP.GammaUp_std_err])
# Q3_data_array.insert_data(tempUpUpError)

# print(Q3_data_array.temp)
# print(Q3_data_array.GammaUp)
# Q3_data_array.plot()