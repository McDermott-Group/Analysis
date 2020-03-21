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


# charge_path = ('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                # '\Parameter\cvg0848kgb_correlation.hdf5')
# P1_path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
            # '{}/General/{}/Charge_resetting/MATLABData/')
            
# CO = ChargeOffset()
# CO.add_dataset(charge_path)
# jumps, sigma = CO.get_jump_sizes()

# date, Q = '03-19-20', 'Q1'
# P1_files = (3,1069+1)
# filenames = [P1_path.format(Q,date) + 'Charge_resetting_{:03d}.mat'.format(i) 
             # for i in range(*P1_files)]
# P1_avg = get_averaged_P1(filenames, '_SB2')
# ax = plot_data(P1_avg, 'File', '{} Avg Excess 1 State'.format(Q))
# ax.plot(np.abs(jumps['Q1']))

# date, Q = '03-19-20', 'Q2'
# P1_files = (3,1069+1)
# filenames = [P1_path.format(Q,date) + 'Charge_resetting_{:03d}.mat'.format(i) 
             # for i in range(*P1_files)]
# P1_avg = get_averaged_P1(filenames, '_SB1')
# ax = plot_data(P1_avg, 'File', '{} Avg Excess 1 State'.format(Q))
# ax.plot(np.abs(jumps['Q2']))

# plt.show()


charge_path = ('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvh0508hpk_correlation.hdf5')
P1_path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
            '{}/General/{}/Charge_resetting/MATLABData/')
            
CO = ChargeOffset()
CO.add_dataset(charge_path)
jumps, sigma = CO.get_jump_sizes()

date, Q = '03-20-20', 'Q1'
P1_files = (4,511+1)
filenames = [P1_path.format(Q,date) + 'Charge_resetting_{:03d}.mat'.format(i) 
             for i in range(*P1_files)]
P1_avg_Q1 = get_averaged_P1(filenames, '_SB2')
ax = plot_data(P1_avg_Q1, 'File', title=Q, label='Avg Excess 1 State')
ax.plot(np.abs(jumps['Q1']), label='Jump Magnitude')
ax.legend()

date, Q = '03-20-20', 'Q2'
P1_files = (4,511+1)
filenames = [P1_path.format(Q,date) + 'Charge_resetting_{:03d}.mat'.format(i) 
             for i in range(*P1_files)]
P1_avg_Q2 = get_averaged_P1(filenames, '_SB1')
ax = plot_data(P1_avg_Q2, 'File', title=Q, label='Avg Excess 1 State')
ax.plot(np.abs(jumps['Q2']), label='Jump Magnitude')
ax.legend()

ax = plot_data(P1_avg_Q1, 'File', 'Avg Excess 1 State', label='Q1')
ax.plot(P1_avg_Q2, label='Q2')

plt.show()