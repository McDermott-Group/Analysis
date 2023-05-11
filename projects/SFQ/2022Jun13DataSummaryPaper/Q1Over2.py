from SFQlib import RB, RB_AllGates, Purity, T1_QP_2D_Linear, T1_QP_2D
from SFQlib import RB_AllGates_Paper
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if 0:  # IRB the optimized result 1.777% error/clifford gate
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    exp_name = 'RB_All'
    date = '2022Jun25Q1Over2V3'  #
    # file_Number = np.arange(0, 6, 1)
    file_Number = np.arange(6, 9, 1)
    RB_file = [
        file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i)
        for i in file_Number]
    # RB_data = RB_AllGates()
    RB_data = RB_AllGates_Paper()
    RB_data.add_data_from_matlab(RB_file)
    RB_data.data_analysis()
    RB_data.plot()

    print('X', RB_data.F_X, RB_data.F_X_std)
    print('Y', RB_data.F_Y, RB_data.F_Y_std)
    print('X/2', RB_data.F_XOver2, RB_data.F_XOver2_std)
    print('-X/2', RB_data.F_XOver2Neg, RB_data.F_XOver2Neg_std)
    print('Y/2', RB_data.F_YOver2, RB_data.F_YOver2_std)
    print('-Y/2', RB_data.F_YOver2Neg, RB_data.F_YOver2Neg_std)

if 1:  # purity data, with 1.137 % incoherence error/clifford
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022Jun22Q1Over2'
    exp_name = 'Purity2D'
    file_Number = np.arange(0, 3, 1)
    file = [
        file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i)
        for i in file_Number]
    Purity_data = Purity()
    Purity_data.add_data_from_matlab(file)
    Purity_data.data_analysis()
    Purity_data.plot()
