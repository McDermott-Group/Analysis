import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
from dataChest import dataChest
from random import randrange


Q_A, Q_B = 2,1
# base_path = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\\'
CO = ChargeOffset()
CO.add_dataset('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvh0508hpk_correlation.hdf5')

path = (r'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\{}\General\03-20-20\Charge_resetting\MATLABData\\')
files = CO.get_files_triggered_on_bad_fit('cvh0508hpk', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
base_file = 'Charge_resetting_00{}.mat'.format(7)

data0 = noiselib.loadmat(path.format('Q'+str(Q_A))+base_file)
SBA0 = np.array(data0['Single_Shot_Occupation_SB{}'.format(Q_A)])
acf = np.zeros(2*50+1)
acf0 = np.zeros(2*50+1)

for f in files:
    data = noiselib.loadmat(f)
    SBA = np.array(data['Single_Shot_Occupation_SB{}'.format(Q_A)])
    SBB = np.array(data['Single_Shot_Occupations_SB{}'.format(Q_B)])
    badIndicies = np.argwhere(np.abs(SBA - SBA0) > 0.2)
    badIndex = badIndicies[0] if len(badIndicies) > 0 else [0]
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(SBA0)
        ax.plot(SBA)
        ax.plot([badIndex],[SBA[badIndex]], 'ro')
        plt.draw()
        plt.pause(0.05)
    r = randrange(10)
    acf += noiselib.crosscorrelate(SBB[badIndex[0]], SBB[badIndex[0]], 50)[0]
    acf0 += noiselib.crosscorrelate(SBB[r], SBB[r], 50)[0]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(acf/len(files), label='Triggered on Q{} charge jump'.format(Q_A))
ax.plot(acf0/len(files), label='Triggered Randomely')
ax.set_xlabel('Lag [100us]')
ax.set_ylabel('Autocorrelation amplitude for Q{}'.format(Q_B))
ax.legend()
plt.draw()
plt.pause(0.05)