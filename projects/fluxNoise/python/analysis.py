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

def plot_data(data, xlabel='', ylabel='', title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.plot(P1_avg)
    plt.draw()
    plt.pause(0.05)
    return ax


date, Q = '03-19-20', 'Q1'
P1_files = (3,1069+1)
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/Charge_resetting/MATLABData/'.format(Q,date))
filenames = [path+'Charge_resetting_{:03d}.mat'.format(i) 
             for i in range(*P1_files)]
P1_avg = get_averaged_P1(filenames, '_SB2')
ax = plot_data(P1_avg, 'File', '{} Avg Excess 1 State'.format(Q))
CO = ChargeOffset()
CO.add_dataset('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvg0848kgb_correlation.hdf5')
jumps, sigma = CO.get_jump_sizes()
ax.plot(np.abs(jumps['Q1']))

date, Q = '03-19-20', 'Q2'
P1_files = (3,1069+1)
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/Charge_resetting/MATLABData/'.format(Q,date))
filenames = [path+'Charge_resetting_{:03d}.mat'.format(i) 
             for i in range(*P1_files)]
P1_avg = get_averaged_P1(filenames, '_SB1')
ax = plot_data(P1_avg, 'File', '{} Avg Excess 1 State'.format(Q))
CO = ChargeOffset()
CO.add_dataset('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvg0848kgb_correlation.hdf5')
jumps, sigma = CO.get_jump_sizes()
ax.plot(np.abs(jumps['Q2']))
plt.show()


raise



# date, Q = '03-17-20', 'Q1'
# P1_files = (0,1003+1)
# QP_files = (0,1003+1)
# filter_levels = (0.05, 0.15) # low, high
# date, Q = '03-17-20', 'Q2'
# P1_files = (0,908+1)
# QP_files = (0,908+1)
# filter_levels = (0.05, 0.05) # low, high
# date, Q = '03-18-20', 'Q3'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.02, 0.02) # low, high
date, Q = '03-18-20', 'Q4'
P1_files = (0,999+1)
QP_files = (0,999+1)
filter_levels = (0.15, 0.2) # low, high

path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/Excess_1_state/MATLABData/'.format(Q,date))
filenames = [path+'Excess_1_state_{:03d}.mat'.format(i) 
             for i in range(*P1_files)]
P1_avg = get_averaged_P1(filenames, '_SB1')
ax = plot_data(P1_avg, 'File', '{} Avg Excess 1 State'.format(Q))

# plot QP tunneling curves
filter_high = P1_avg > filter_levels[1]
filter_low = P1_avg < filter_levels[1]
ax.plot([0, P1_avg.size], [filter_levels[0]]*2)
ax.plot([0, P1_avg.size], [filter_levels[1]]*2)
file_indicies = np.arange(*QP_files)
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/QP_Tunneling_PSD/MATLABData/'.format(Q,date))
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in file_indicies[filter_high]]
QPT = QPTunneling()
QPT.add_datasets(filenames)
figN = QPT.plot_psd(label=Q+' high')
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in file_indicies[filter_low]]
QPT = QPTunneling()
QPT.add_datasets(filenames)
QPT.plot_psd(figNum=figN, label=Q+' low')

# path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        # 'Q1/General/03-17-20/QP_Tunneling_PSD/MATLABData/')
# QPT = QPTunneling()
# filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) for i in range(100)]
# QPT.add_datasets(filenames)
# QPT.plot_psd()