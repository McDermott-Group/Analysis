import numpy as np
import matplotlib.pyplot as plt
import noiselib
from QPTunneling import *
from ChargeOffset import *

def apply_infidelity_correction(o, n_bins=7, thresh=0.5):
    for trial in o:
        # trial[...] = movingmean(trial, n_bins) > thresh
        a = trial.astype(np.float)
        for i in range(len(a)):
            i_0 = max(0,i-int(np.floor(n_bins/2.)))
            i_end = min(len(a), i+int(np.ceil(n_bins/2.)))
            a[i] = np.mean(a[ i_0:i_end ])
        trial[...] = a > thresh
        
def apply_infidelity_correction2(o, n_bins=7, thresh=0.5):
    a = o.astype(np.float)
    for i in range(len(a)):
        i_0 = max(0,i-int(np.floor(n_bins/2.)))
        i_end = min(len(a), i+int(np.ceil(n_bins/2.)))
        a[i] = np.mean(a[ i_0:i_end ])
    o[:] = a > thresh

def get_averaged_P1(files, label=''):
    P1_avg = np.array([])
    for f in files:
        data = noiselib.loadmat(f)
        o = np.array(data['Single_Shot_Occupation{}'.format(label)])
        P1_avg = np.append(P1_avg, o)
    return P1_avg
    
def get_averaged_flips(files, label=''):
    flip_avg = np.array([])
    for f in files:
        data = noiselib.loadmat(f)
        o = np.array(data['Single_Shot_Occupations{}'.format(label)])
        apply_infidelity_correction(o, 9)
        flip_avg = np.append(flip_avg, np.mean(np.abs(np.diff(o, axis=1))))
    return flip_avg

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


###fluxNoise
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
# date, Q = '03-18-20', 'Q4'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.15, 0.2) # low, high
### fluxNoise2
# date, Q = '03-17-20', 'Q1'
# P1_files = (0,1003+1)
# QP_files = (0,1003+1)
# filter_levels = (0.05, 0.15) # low, high
date, Q = '04-10-20', 'Q2'
P1_files = (0,999+1)
QP_files = (0,999+1)
filter_levels = (0.05, 0.05) # low, high

# plot P1
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/Excess_1_state/MATLABData/'.format(Q,date))
filenames = [path+'Excess_1_state_{:03d}.mat'.format(i)
             for i in range(*P1_files)]
P1_avg = get_averaged_P1(filenames, '_SB1')
ax = plot_data(noiselib.movingmean(P1_avg, 50), 'File', '{} Avg Excess 1 State'.format(Q))

file_indicies = np.arange(*QP_files)
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/QP_Tunneling_PSD/MATLABData/'.format(Q,date))
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in file_indicies]
flip_avg = get_averaged_flips(filenames, '_SB1')
ax.plot(noiselib.movingmean(flip_avg, 50), label='Avg Flip Rate')
plt.draw()
plt.pause(0.05)