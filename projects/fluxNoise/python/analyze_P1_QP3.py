import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from QPTunneling import *
from ChargeOffset import *

def get_averaged_P1(files, label=''):
    P1_avg = np.array([])
    for f in files:
        try:
            data = noiselib.loadmat(f)
        except:
            print 'corrupted:', f
            data = { 'Single_Shot_Occupation{}'.format(label): np.nan }
        o = np.array(data['Single_Shot_Occupation{}'.format(label)])
        P1_avg = np.append(P1_avg, o)
    return P1_avg
    
def get_averaged_flips(files, label=''):
    flip_avg = np.array([])
    for f in files:
        data = noiselib.loadmat(f)
        o = np.array(data['Single_Shot_Occupations{}'.format(label)])
        o = noiselib.apply_infidelity_correction_HMM(o, fidelity=[0.95,0.75])
        # o = noiselib.apply_infidelity_correction(o, 9)
        flip_avg = np.append(flip_avg, np.mean(np.abs(np.diff(o, axis=1))))
    return flip_avg

def plot_data(data, xlabel='', ylabel='', title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.plot(data)
    plt.draw()
    plt.pause(0.05)
    return ax


###fluxNoise2
# date, Q = '04-28-20', 'Q3'
# P1_files = (0,821+1)
# QP_files = (0,821+1)
# date, Q = '04-29-20', 'Q3'
# P1_files = (0,110+1)
# QP_files = (0,110+1)
# date, Q = '04-30-20', 'Q4'
# P1_files = (234,470+1)
# QP_files = (1,237+1)
date, Q = '05-01-20', 'Q3'
# P1_files = (0,664+1)
# QP_files = (0,664+1)
P1_files = (0,6+1)
QP_files = (0,6+1)

# plot P1
path = ('Z:/mcdermott-group/data/fluxNoise2/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/{}/MATLABData/'.format(Q,date,'{}'))
filenames_I = [path.format('P1_I') + 'P1_I_{:03d}.mat'.format(i)
                for i in range(*P1_files)]
filenames_X = [path.format('P1_X') + 'P1_X_{:03d}.mat'.format(i)
                for i in range(*P1_files)]
P1_avg_I = get_averaged_P1(filenames_I, '')
P1_avg_X = get_averaged_P1(filenames_X, '')
ax = plot_data(noiselib.movingmean(np.transpose([P1_avg_I, P1_avg_X]), 1), 'File', '{} Avg Excess 1 State'.format(Q))

path = ('Z:/mcdermott-group/data/fluxNoise2/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/QP_Tunneling_PSD/MATLABData/'.format(Q,date))
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in np.arange(*QP_files)]
flip_avg = get_averaged_flips(filenames, '')
ax.plot(20*noiselib.movingmean(flip_avg, 50), label='20 x Avg Flip Rate')
plt.draw()
# plt.pause(0.05)
plt.show()