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
# date = '03-19-20'
# P1_files = {'Q1':(3,1069+1), 'Q2':(3,1069+1)}


charge_path = ('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvh0508hpk_correlation.hdf5')
P1_path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
            '{}/General/{}/Charge_resetting/MATLABData/')
date = '03-20-20'
P1_files = {'Q1':(4,511+1), 'Q2':(4,511+1)}
            
            
CO = ChargeOffset()
CO.add_dataset(charge_path)
jumps, sigma = CO.get_jump_sizes()
P1_avg = {}
for QA,QB in [('Q1','Q2'), ('Q2','Q1')]:
    filenames = [P1_path.format(QA,date) + 'Charge_resetting_{:03d}.mat'.format(i) 
                 for i in range(*P1_files[QA])]
    P1_avg[QB] = get_averaged_P1(filenames, '_SB{}'.format(QB[-1]))
    ax = plot_data(P1_avg[QB], 'File', title=QB, label='Avg Excess 1 State')
    ax.plot(np.abs(jumps[QB]), label='Jump Magnitude')
    ax.legend()
    plt.draw()
    plt.pause(0.05)

ax = plot_data(P1_avg['Q1'], 'File', 'Avg Excess 1 State', label='Q1')
ax.plot(P1_avg['Q2'], label='Q2')
plt.draw()
plt.pause(0.05)