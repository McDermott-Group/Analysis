"""
Analyze the up transition rate based on SFQ pulses
"""

from SFQlib import QP_Up
import numpy as np


# """
# P1 1D Fit
# """
file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
date = '01-10-22'
# experiment_name = 'SFQ_OffResonantDrive_Up_Study_1us'
experiment_name = 'SFQ_Rabi_3GHz_Drive'

file_Number = [0]
P1_file = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_data = QP_Up()
P1_data.add_data_from_matlab(P1_file)
P1_data.plot()


