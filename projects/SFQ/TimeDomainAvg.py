from SFQlib import DeepSubharmonics
import matplotlib.pyplot as plt
import numpy as np

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')

# date = '2022Jun09_PoisonTimeDomain'
# date = '06-22-22'
# date = '2022Jun11_PoisonTimeDomainQ1'
# date = '2022Jun21SFQ1PoisonQ1Q2P1'
date = '2022Jun22SFQ1PoisonQ1Q2P1'
file_Number = np.arange(0, 10, 1)
# file_Number = np.arange(0, 9, 1)
# file_Number = np.arange(0, 50, 1)
# file_Number = np.arange(11, 16, 1)
# file_Number = [0, 2, 4, 5]
# files have issue: 4
# experiment_name = ('Rabi_SFQ_PoisonIdle_First50ns')
# experiment_name = ('Rabi_SFQ_NoPoison_7dB')
# experiment_name = ('Rabi_SFQ_Poison_7dB_Logns')
# experiment_name = ('Rabi_SFQ_Poison_7dB_ShortLong')
# experiment_name = ('Rabi_SFQ_Poison')
# experiment_name = ('Rabi_SFQ_Poison_IBias2Q1')
# experiment_name = ('SFQ1Poison_Q1P1_V1Long')
# experiment_name = ('SFQ1Poison_Q2P1_V1Long')
# experiment_name = ('SFQ1Poison_Q2P1_V1Short')
# experiment_name = ('SFQ1Poison7GHzQ1_Test')
# experiment_name = ('SFQ1Poison7GHzQ1_Test_13dBm')
# experiment_name = ('SFQ1Poison7GHzQ2_Test')
# experiment_name = ('PoisonP1LongRO')
# experiment_name = ('PoisonP1Q1_Test')
# experiment_name = ('PoisonP1Q2_V1')
# experiment_name = ('PoisonP1Q1_V1')
experiment_name = ('PoisonP1Q2_V1_LongCorrect')
# experiment_name = ('PoisonP1Q1_V1_Long')
file_list_Q1 = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
P1 = DeepSubharmonics()
P1.add_data_from_matlab(file_list_Q1)
P1.plot()