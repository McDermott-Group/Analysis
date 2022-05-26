from SFQlib import RB, RB_AllGates
import numpy as np


# file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# # date = '2022May20GateTune'
# date = '2022May21GateTune'
# # exp_name = 'RB_AfterCalStats'
# exp_name = 'RB'
# n = 0
# file_Number = [0+n, 6+n, 12+n, 18+n, 24+n, 30+n, 36+n]
# RB_file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
# RB_data = RB()
# RB_data.add_data_from_matlab(RB_file)
# RB_data.data_analysis()
# RB_data.plot()

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
date = '05-26-22'
# exp_name = 'RB_AfterCalStats'
exp_name = 'RB_OneRef'
file_Number = [11, 12, 13]
RB_file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
RB_data = RB_AllGates()
RB_data.add_data_from_matlab(RB_file)
RB_data.data_analysis()
RB_data.plot()