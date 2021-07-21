import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from QPTunneling import *
import ChargeOffset
importlib.reload(ChargeOffset)
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


# base_path = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'cvi1713bxq_correlation.hdf5')
# # CO.plot_charge_offset()
# # CO.plot_jump_sizes()
# path = r'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\03-21-20\Charge_resetting_QP\HDF5Data'
# base_file = 'cvi1721ckw_Charge_resetting_QP.hdf5'
# files = CO.get_files_triggered_on_bad_fit('cvi1713bxq', 'Q1', path)
# var_list = ['Single Shot Occupation SB{}'.format(i+1) for i in range(2)]

# dc = dataChest(path.split('\\'))
# dc.openDataset(base_file)
# data0 = dc.getData(
    # variablesList=var_list)
# data0 = data0.transpose()

# for f in files[:3]:
    # dc.openDataset(f)
    # data = dc.getData(
        # variablesList=var_list)
    # data = data.transpose()
    # badIndicies = np.argwhere(np.abs(data[0] - data0[0]) > 0.2)
    # badIndex = badIndicies[0] if len(badIndicies) > 0 else 0
    # # fig = plt.figure()
    # # ax = fig.add_subplot(111)
    # # ax.plot(data0[0])
    # # ax.plot(data[0])
    # # ax.plot([badIndex],[data[0][badIndex]], 'ro')
    # # plt.draw()
    # # plt.pause(0.05)
    # data = dc.getData(
        # variablesList=['Single Shot Occupations SB2'])
    # data = data.transpose()
    # print data.shape
# plt.show()