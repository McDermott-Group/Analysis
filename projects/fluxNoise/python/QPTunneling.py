import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt

class QPTunneling(object):
    
    def __init__(self, fs=1/100e-6):
        self.fs = fs
        vis = 1
        self.transfer = 1/vis**2
    
    def add_datasets(self, file_path):
        if type(file_path) == str:
            file_path = [file_path]
        for f in file_path:
            data = noiselib.loadmat(f)
            o = np.array(data['Single_Shot_Occupations_SB1'])
            trials, reps = o.shape
            if not hasattr(self, 'data'):
                self.data = np.empty((0,reps))
            self.data = np.concatenate([self.data, o])
    
    def get_psd(self, window_averaging = True):
        cpsd, f = noiselib.partition_and_avg_psd(self.data, self.fs)
        cpsd = self.transfer * cpsd
        if window_averaging:
            cpsd = noiselib.window_averaging(cpsd)
        return cpsd, f
        
    def plot_psd(self, window_averaging = True, figNum=None, label=None):
        psd, f = self.get_psd(window_averaging)
        fig = plt.figure(figNum)
        ax = fig.add_subplot(111)
        ax.set_title('PSD')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('$S_\eta (\eta^2/Hz)$')
        ax.plot(f, psd, label=label)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        ax.grid()
        plt.draw()
        plt.pause(0.05)
        return fig.number

# Q = 'Q1'
# date = '03-17-20'
# P1_files = (0,1003+1)
# QP_files = (0,1003+1)
# filter_levels = (0.05, 0.15) # low, high
# Q = 'Q2'
# date = '03-17-20'
# P1_files = (0,908+1)
# QP_files = (0,908+1)
# filter_levels = (0.05, 0.05) # low, high
# Q = 'Q3'
# date = '03-18-20'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.02, 0.02) # low, high
Q = 'Q4'
date = '03-18-20'
P1_files = (0,999+1)
QP_files = (0,999+1)
filter_levels = (0.17, 0.17) # low, high
path = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/'
        '{}/General/{}/Excess_1_state/MATLABData/'.format(Q,date))
P1_avg = []
for i in range(*P1_files):
    data = noiselib.loadmat(path+'Excess_1_state_{:03d}.mat'.format(i))
    o = np.array(data['Single_Shot_Occupations_SB1'])
    P1_avg.append(np.mean(o))
P1_avg = np.array(P1_avg)
# plot P1
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('')
ax.set_xlabel('File')
ax.set_ylabel('{} Avg Excess 1 State'.format(Q))
ax.plot(P1_avg)
plt.draw()
plt.pause(0.05)
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