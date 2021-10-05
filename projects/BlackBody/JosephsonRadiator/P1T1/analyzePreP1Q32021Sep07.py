from antennalib import GMMFit
import numpy as np

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = '09-11-21'

experiment = ('Q3_P1_17500uDACJ2')

IName = experiment+'_I_Cal'
XName = experiment+'_X_Cal'
J2Name = experiment

file_Number = np.arange(0, 1, 1)
I_file_list = [file_path.format(date, IName, IName) + '_{:03d}.mat'.format(i) for i in file_Number]
X_file_list = [file_path.format(date, XName, XName) + '_{:03d}.mat'.format(i) for i in file_Number]
J2_file_list = [file_path.format(date, J2Name, J2Name) + '_{:03d}.mat'.format(i) for i in file_Number]

Q3 = GMMFit()
Q3.add_data_and_process_data(path=[I_file_list, X_file_list, J2_file_list])
Q3.plot_data(s='I')
Q3.plot_data(s='X')
Q3.plot_data(s='J', show=True)