from SFQlib import DeepSubharmonics
import matplotlib.pyplot as plt
import numpy as np

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')

# date = '2022Apr26SFQStep'
# date = '2022Apr27SFQStep'
# date = '2022Apr27SFQStepReal'
# date = '2022May09sSFQStepReal'
# date = '2022May10SFQStepReal'
# date = '2022May11SFQStepReal_v1'
date = '2022May12SFQStepReal'
# date = '04-27-22'
# Z:\mcdermott-group\data\sfq\MCM_NIST\LIU\MCM13\2022May09sSFQStepReal\SFQ_SB_1D_Steps

# file_Number = [0, 1, 3, 4, 6, 12, 13, 14, 15, 16, 17, 18, 19, 23, 24, 26, 27,
#                28, 30, 32, 33, 34, 35, 37, 39, 41, 43, 45, 49]
file_Number = np.arange(0, 14, 1)
# file_Number = np.arange(25, 50, 1)
experiment_name = ('SFQ_SB_1D_Steps')
file_list_Q1 = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
P1 = DeepSubharmonics()
P1.add_data_from_matlab(file_list_Q1)
P1.plot()