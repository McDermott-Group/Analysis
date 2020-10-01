import numpy as np
import matplotlib.pyplot as plt
import noiselib
from QPTunneling import *
from ChargeOffset import *

def get_averaged_P1(files, label=''):
    P1_avg = np.array([])
    for f in files:
        data = noiselib.loadmat(f)
        o = np.array(data['Single_Shot_Occupation{}'.format(label)])
        P1_avg = np.append(P1_avg, np.mean(o))
    return P1_avg

def plot_data(data, xlabel='', ylabel='', title='', label=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.plot(data, label=label)
    plt.draw()
    plt.pause(0.05)
    return ax


charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvi1713bxq_correlation.hdf5')
QP_path = ('Z:/mcdermott-group/data/fluxNoise2/DR1 - 2019-12-17/CorrFar/'
            'Q1Q2Corr/General/{}/Charge_resetting_QP/MATLABData/')
            
CO = ChargeOffset()
CO.add_dataset(charge_path)
jumps, sigma = CO.get_jump_sizes()
large_jumps = {l: np.append(False, np.abs(jumps[l]) > 0.1) for l in jumps.keys()}

date = '03-21-20'
QP_files_Q1 = np.arange(0, 489+1, 2)
QP_files_Q2 = np.arange(1, 489+1, 2)
# QP_files_Q1 = np.arange(0, 49+1, 2)
# QP_files_Q2 = np.arange(1, 49+1, 2)
filenames_Q1_jump = [QP_path.format(date) + 'Charge_resetting_QP_{:03d}.mat'.format(i) 
                    for i in QP_files_Q1[large_jumps['Q1']]]
filenames_Q1_nojump = [QP_path.format(date) + 'Charge_resetting_QP_{:03d}.mat'.format(i) 
                        for i in QP_files_Q1[~large_jumps['Q1']]]
filenames_Q2_jump = [QP_path.format(date) + 'Charge_resetting_QP_{:03d}.mat'.format(i) 
                    for i in QP_files_Q2[large_jumps['Q2']]]
filenames_Q2_nojump = [QP_path.format(date) + 'Charge_resetting_QP_{:03d}.mat'.format(i) 
                        for i in QP_files_Q2[~large_jumps['Q2']]]

# plot QP tunneling curves
for i,file_list in enumerate([filenames_Q1_jump, filenames_Q1_nojump, filenames_Q2_jump, filenames_Q2_nojump]):
    QPT = QPTunneling()
    QPT.add_datasets(file_list, 'Single_Shot_Occupations_SB{}'.format([2,2,1,1][i]))
    QPT.plot_psd(figNum=111)#, label=Q+' high')

plt.show()