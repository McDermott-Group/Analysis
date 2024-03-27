from .SFQlib import RB, RB_AllGates, Purity
import numpy as np


file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
# date = '06-03-22'
# date = '2022Jun03_Purity'
date = '2022Jun05_PurityVsRB'
# exp_name = 'RB_Purity_New'
# exp_name = 'RB_Purity'
# exp_name = 'PurityNorm'
exp_name = 'Purity2D'
n = 7
# file_Number = [0]
file_Number = np.arange(0, 10, 1)
# file_Number = [2, 3, 4, 6, 7]
# file_Number = [1, 2, 3, 4, 5, 6, 7, 8]
# file_Number = [6+n, 12+n, 18+n, 24+n, 30+n]
# file_Number = [0+n, 6+n, 12+n]
file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
Purity_data = Purity()
Purity_data.add_data_from_matlab(file)
Purity_data.data_analysis()
# print('n=', n)
# print('Fidelity and std=[ {:.3f}, {:.3f} ]%'.format(RB_data.F*100, RB_data.F_std*100))
Purity_data.plot()
