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


Q_A, Q_B = 1,2
CO = ChargeOffset()
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\cvi1713bxq_correlation.hdf5')
path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\03-21-20\Charge_resetting_QP\MATLABData\\'
files = CO.get_files_triggered_on_bad_fit('cvi1713bxq', 'Q{}'.format(Q_A), path)
for i,f in enumerate(files):
    beginning, end = f.rsplit('_',1)
    num, ext = end.split('.',1)
    files[i] = beginning + '_{:03d}.'.format(int(num)+Q_A-1) + ext

base_file = 'Charge_resetting_QP_00{}.mat'.format(7+Q_A)
data0 = noiselib.loadmat(path+base_file)
SBA0 = np.array(data0['Single_Shot_Occupation_SB{}'.format(Q_A)])
QPT0 = QPTunneling()
QPT = QPTunneling()
acf = np.zeros(2*50+1)
acf0 = np.zeros(2*50+1)

for f in files:
    print('f=', f)  # V
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
    # print(SBB[badIndex]) #V
    QPT.add_data(SBB[badIndex])
    r = randrange(10)
    QPT0.add_data(SBB[[r]])
    acf += noiselib.crosscorrelate(SBB[badIndex[0]], SBB[badIndex[0]], 50)[0]
    acf0 += noiselib.crosscorrelate(SBB[r], SBB[r], 50)[0]
QPT.plot_psd(figNum=111)
QPT0.plot_psd(figNum=111)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(acf)
ax.plot(acf0)
plt.draw()
plt.pause(0.05)
plt.show()