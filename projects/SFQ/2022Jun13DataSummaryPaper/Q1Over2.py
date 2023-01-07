from SFQlib import RB, RB_AllGates, Purity, T1_QP_2D_Linear, T1_QP_2D
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if 0:  # IRB the optimized result 1.777% error/clifford gate
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022Jun25Q1Over2V3'  #
    exp_name = 'RB_All'
    file_Number = np.arange(6, 9, 1)
    RB_file = [
        file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i)
        for i in file_Number]
    RB_data = RB_AllGates()
    RB_data.add_data_from_matlab(RB_file)
    RB_data.data_analysis()
    RB_data.plot()

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
