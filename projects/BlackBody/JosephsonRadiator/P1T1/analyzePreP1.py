from antennalib import QP_Up, Up_array
import numpy as np

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorP1Pre_2021Sep01'

experiment_name = ('Q1_P1_0uDACJ2')
# experiment_name = ('Q4_P1_54300uDACJ2')

file_Number_Tot = np.arange(0, 20, 1)
file_Number_Extract = np.array([])
file_Number = np.setdiff1d(file_Number_Tot, file_Number_Extract)
file_list = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]

QP = QP_Up()
QP.add_data_from_matlab(file_list, data_type1='Occupation_Filtered')
QP.plot_Chi()
QP.plot()