from SFQlib import DeepSubharmonics
import matplotlib.pyplot as plt
import numpy as np

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')


date = '2022Jun28PoisonTimeDomainRealV1'
# file_Number = np.arange(0, 300, 1)
file_Number = np.arange(0, 150, 1)
# experiment_name = ('PoisonQ1')
experiment_name = ('PoisonQ2')
file_list_Q1 = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
P1 = DeepSubharmonics()
P1.add_data_from_matlab(file_list_Q1)
P1.plot()