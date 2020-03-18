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
        
    def plot_psd(self, window_averaging = True, figNum=None):
        psd, f = self.get_psd(window_averaging)
        fig = plt.figure(figNum)
        ax = fig.add_subplot(111)
        ax.set_title('PSD')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('$S_\eta (\eta^2/Hz)$')
        ax.plot(f, psd)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.grid()
        plt.draw()
        plt.pause(0.05)
        return fig.number

path = 'Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q2/General/03-17-20/Excess_1_state/MATLABData/'
P1_avg = []
for i in range(908):
    data = noiselib.loadmat(path+'Excess_1_state_{:03d}.mat'.format(i))
    o = np.array(data['Single_Shot_Occupations_SB1'])
    P1_avg.append(np.mean(o))
# plot P1
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('')
ax.set_xlabel('File')
ax.set_ylabel('Avg Excess 1 State')
ax.plot(P1_avg)
plt.draw()
plt.pause(0.05)
# plot QP tunneling curves
P1_avg = np.array(P1_avg)
filter_high = P1_avg > 0.15
filter_low = P1_avg < 0.05
file_indicies = np.arange(P1_avg.size)
path = 'Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q2/General/03-17-20/QP_Tunneling_PSD/MATLABData/'
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in file_indicies[filter_high]]
QPT = QPTunneling()
QPT.add_datasets(filenames)
figN = QPT.plot_psd()
filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) 
            for i in file_indicies[filter_low]]
QPT = QPTunneling()
QPT.add_datasets(filenames)
QPT.plot_psd(figNum=figN)

# path = 'Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q1/General/03-17-20/QP_Tunneling_PSD/MATLABData/'
# QPT = QPTunneling()
# filenames = [path+'QP_Tunneling_PSD_{:03d}.mat'.format(i) for i in range(100)]
# QPT.add_datasets(filenames)
# QPT.plot_psd()