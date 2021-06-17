from antennalib import QP_Up, Up_array
import numpy as np

Q3_data_array = Up_array()

QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')


date = 'P1Pre2021Jun11'
experiment_name_global = 'Q3_P1_'

# temp_list = [78, 200, 426]
# temp_list = [78, 185, 194, 200, 237, 274, 306, 336, 365, 397, 426]
temp_list = [78, 194, 200, 237, 274, 306, 336, 365, 397, 426]
for temp in temp_list:
    experiment_name = experiment_name_global + str(temp) + 'mk'

    file_Number_Tot = np.arange(0, 50, 1)
    file_Number_Extract = np.array([])
    file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
    #
    QP_file = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
    QP = QP_Up()
    QP.add_data_from_matlab(QP_file, temp)
    # QP.plot()
    tempUpUpError = np.array([temp, QP.GammaUp, QP.GammaUp_std_err])
    Q3_data_array.insert_data(tempUpUpError)
    # print(temp, QP.GammaUp, QP.GammaUp_std_err)
    # print("[%3d, %4d,%4d]"%(temp, QP.GammaUp, QP.GammaUp_std_err))

# print(Q3_data_array.temp)
# print(Q3_data_array.GammaUp)
Q3_data_array.plot()