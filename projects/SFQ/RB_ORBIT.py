from SFQlib import RB, RB_AllGates, Purity
import numpy as np


# file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# # date = '2022May20GateTune'
# # date = '2022May21GateTune'
# # date = '2022May29RBCheck'
# date = '2022May30RB'
# # exp_name = 'RB_AfterCalStats'
# exp_name = 'RB_Coherence'
# n = 7
# file_Number = [0+n, 8+n, 16+n, 24+n, 32+n, 40+n, 48+n, 56+n, 64+n, 72+n, 80+n]
# # file_Number = [6+n, 12+n, 18+n, 24+n, 30+n]
# # file_Number = [0+n, 6+n, 12+n]
# RB_file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
# RB_data = RB()
# RB_data.add_data_from_matlab(RB_file)
# RB_data.data_analysis()
# print('n=', n)
# print('Fidelity and std=[ {:.3f}, {:.3f} ]%'.format(RB_data.F*100, RB_data.F_std*100))
# RB_data.plot()

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# date = '06-23-22'
# date = '2022Jun17Q2'
# date = '06-19-22'
# date = '2022May30RB'    # the optimized result 1.19% error/clifford gate
# date = '2022Jun03_Purity'
# date = '2022Jun05_PurityVsRB'
# date = '2022Jun19Q2'
# date = '2022Jun20Over39nQ1'
date = '2022Jun23Over2'
# exp_name = 'RB_AfterCalStats'
exp_name = 'RB_All'
# exp_name = 'Purity2D'
# file_Number = [0, 1, 2]
file_Number = np.arange(5, 10, 1)
RB_file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
RB_data = RB_AllGates()
RB_data.add_data_from_matlab(RB_file)
RB_data.data_analysis()
RB_data.plot()

# Purity_data = Purity()
# Purity_data.add_data_from_matlab(RB_file)
# Purity_data.data_analysis()
# Purity_data.plot()