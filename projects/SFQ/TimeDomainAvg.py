from SFQlib import DeepSubharmonics
import matplotlib.pyplot as plt
import numpy as np

file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')

# date = '06-26-22'
date = '06-27-22'
# date = '2022Jun26PoisonTimeDomain'
# date = '2022Jun27PoisonTimeDomain'
file_Number = np.arange(0, 25, 1)
# experiment_name = ('PoisonP1ShortQ1')
# experiment_name = ('PoisonP1ShortFineQ1')
# experiment_name = ('PoisonP1ShortFineQ2')
# experiment_name = ('PoisonP1ShortFine96nsQ2')
# experiment_name = ('PoisonP1ShortFine96nsQ1')
# experiment_name = ('PoisonP1ShortFine96nsQ1_50nsPoison')
# experiment_name = ('PoisonP1ShortFine96nsQ2_50nsPoison')
# experiment_name = ('Poison50nsRO96nsQ1')
experiment_name = ('Poison50nsRO96nsQ2')
file_list_Q1 = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
P1 = DeepSubharmonics()
P1.add_data_from_matlab(file_list_Q1)
P1.plot()