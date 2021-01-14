import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class QPTunneling(object):
    
    def __init__(self, fs=1/100e-6):
        self.fs = fs
        vis = 1
        self.transfer = 1/vis**2
        self.n_rows = 0
    
    def add_datasets(self, file_path, data_str='Single_Shot_Occupations'):
        if type(file_path) == str:
            file_path = [file_path]
        for f in file_path:
            data = noiselib.loadmat(f)
            o = np.array(data[data_str])
            # o = noiselib.apply_infidelity_correction(o, 9)
            self.add_data(o)
    
    def add_data(self, o):
        self.n_rows += 1 # don't divide by number of rows because these are already averaged in partition_and_avg_psd
        cpsd, f = noiselib.partition_and_avg_psd(o, self.fs)
        if not hasattr(self, 'cpsd'):
            self.cpsd, self.f = cpsd, f
        else:
            self.cpsd += cpsd
    
    def get_psd(self, window_averaging = True):
        cpsd = self.transfer * self.cpsd / self.n_rows
        if window_averaging:
            cpsd = noiselib.window_averaging(cpsd)
        return cpsd, self.f
    
    def get_fit(self, window_averaging = True):
        def y(x, a, gamma):
            return gamma*a/(np.pi**2*x**2 + gamma**2)
        psd, f = self.get_psd(window_averaging)
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]
        popt, pcov = curve_fit(y, f, psd, bounds=(0,np.inf))
        print popt
        return popt[-1], y(f, *popt), f
        
    def plot_psd(self, window_averaging = True, figNum=None, label=None, fit=True):
        fig = plt.figure(figNum)
        ax = fig.add_subplot(111)
        ax.set_title('PSD')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('$S_\eta (\eta^2/Hz)$')
        
        if label is None:
            label = 'psd'
        psd, f = self.get_psd(window_averaging)
        ax.plot(f, psd, label=label)
        if fit:
            rate, psd_fit, f_fit = self.get_fit(window_averaging)
            ax.plot(f_fit, psd_fit, '--', label='{} fit [{:.2f}Hz]'.format(label,rate))
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        ax.grid()
        plt.draw()
        plt.pause(0.05)
        return fig.number