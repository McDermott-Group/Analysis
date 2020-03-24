import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit

class ChargeOffset(object):

    def __init__(self):
        self.file_tags = []
        self.labels = {}
        self.time = {}
        self.unwrapped_charge = {}
    
    def analyze_charge_jump_correlations(self, file_path):
        pass
    
    def add_dataset(self, file_path):
        fpath, fname = os.path.split(file_path)
        ftag = fname.split('_')[0]
        self.file_tags += [ftag]
        dc = dataChest(fpath.split('\\'))
        dc.openDataset(fname)
        
        varsList = dc.getVariables()
        timeVar = varsList[0][0][0]
        offsetVars = [var[0] for var in varsList[1] if 'R2' not in var[0]]
        data = dc.getData( variablesList=[timeVar]+offsetVars )
        data = data.transpose()
        labels = [l[-2:] for l in offsetVars]
        
        time = data[0][1:] - data[0][0]
        period2e = {l: dc.getParameter('2e Period {}'.format(l)) for l in labels}
        unwrapped_charge = { l: noiselib.unwrap_voltage_to_charge(
                                data[i+1], period2e[l]/2, 2/period2e[l] )
                                for i,l in enumerate(labels) }
        unwrapped_charge = { l: unwrapped_charge[l] - unwrapped_charge[l][0] 
                                for l in labels } # zero the trace
        # charge_jump_sizes = {l: unwrapped_charge[l][1:] - unwrapped_charge[l][:-1]
                                # for l in labels}
        # large_jump[ftag]s = {l:charge_jump_sizes[l] > 0.1 for l in labels}
        
        self.labels[ftag] = labels
        self.time[ftag] = time
        self.unwrapped_charge[ftag] = unwrapped_charge
        self._last_dataset_added = ftag
        return self
    
    def limit_dataset(self, label, dataset=None, start=None, end=None):
        """limits the data for label by setting the values in [start,end) to nan"""
        if dataset is None:
            dataset = self._last_dataset_added
        self.unwrapped_charge[dataset][label][start:end] = np.nan
        return self
        
    def _above(self, array, threshold=0.1):
        """Equivilent to array > threshold, but deals with nans"""
        array = array.copy()
        jump_size = np.abs(array, array, where=~np.isnan(array))
        is_above = np.greater( jump_size, threshold, jump_size, where=~np.isnan(array) )
        is_above[np.isnan(is_above)] = False
        return is_above.astype(bool)
    
    def get_charge_offset(self, datasets = None):
        if datasets is None:
            datasets = self.file_tags
        time = np.array([0])
        offset = {}
        for ftag in datasets:
            n_old = time.size
            time = np.append( time, time[-1] + self.time[ftag] )
            n_new = self.time[ftag].size
            for l in self.labels[ftag]:
                # fill with nan for length of other traces
                new_charge = self.unwrapped_charge[ftag][l]
                if l not in offset:
                    offset[l] = np.append( np.full(n_old-1, np.nan), [0] )
                lastValIndex = np.argwhere(~np.isnan(offset[l]))[-1][0]
                offset[l] = np.append( offset[l], offset[l][lastValIndex] + new_charge )
            for l in set(offset.keys()) - set(self.labels[ftag]):
                offset[ftag] = np.append( offset[ftag], np.full(n_new, np.nan) )
        return time, offset
    
    def get_jump_sizes(self, datasets = None):
        time, offset = self.get_charge_offset(datasets)
        jumps = {l: offset[l][1:] - offset[l][:-1] for l in offset.keys()}
        
        # find 2*sigma jump size for each label
        def gaus(x, b, A, mu, sigma):
            return b + A * exp(-(x-mu)**2/(2*sigma**2))
        sigma = {}
        for l in jumps.keys():
            h, bins = np.histogram(jumps[l][~np.isnan(jumps[l])], bins=300)
            x = (bins[1:]+bins[:-1])/2
            popt, pcov = curve_fit(gaus, x, h, p0=[0,h.max(),0,0.1])
            sigma[l] = popt[-1]
        
        return jumps, sigma
        
    def plot_charge_offset(self, datasets = None):
        time, offset = self.get_charge_offset(datasets)
        jumps, sigma = self.get_jump_sizes(datasets)
        largeJumps = {l: np.append(self._above(jumps[l], 2*sigma[l]), True) | 
                         np.append(True, self._above(jumps[l], 2*sigma[l]))
                        for l in jumps.keys()}
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Charge Offset')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Charge Offset [e]')
        for l in offset.keys():
            ax.plot(time, offset[l], label=l)
            ax.plot(time, np.ma.masked_where(~largeJumps[l], offset[l]), label=l)
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        
    def plot_jump_sizes(self, datasets = None):
        jumps, sigma = self.get_jump_sizes(datasets)
        
        for l in jumps.keys():
            n_tot = np.sum(~np.isnan(jumps[l]))
            n_big = np.sum(self._above(jumps[l], 2*sigma[l]))
            print('{}: {} / {} = {:.2f}% = {:.2f}*T between jumps'.format(
                    l, n_big, n_tot, 100.*n_big/n_tot, 1.*n_tot/n_big))
        print( {l: '{:.3f}'.format(2*s) for l,s in sigma.items()} )
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Jump Sizes')
        ax.set_xlabel('N')
        ax.set_ylabel('Jump Size [e]')
        for l in jumps.keys():
            ax.plot(np.abs(jumps[l], jumps[l], where=~np.isnan(jumps[l])), label=l)
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        
    def plot_charge_correlation(self, label1, label2, datasets=None):
        jumps, sigma = self.get_jump_sizes(datasets)
        jumps1 = jumps[label1]
        jumps2 = jumps[label2]
        corrJumps = self._above(jumps1, 2*sigma[label1]) & self._above(jumps2, 2*sigma[label2])
        
        print( 'Correlated Jumps: {}-{}'.format(label1,label2) )
        for l,j in [(label1,jumps1),(label2,jumps2)]:
            print( '    {} / {} = {:.2f}%'.format( 
                    np.sum(corrJumps), 
                    np.sum(self._above(j,2*sigma[l])),
                    100.*np.sum(corrJumps)/np.sum(self._above(j,2*sigma[l])) ) )
                
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,0.5)
        ax.set_title('Charge Jumps')
        ax.set_xlabel('Jump Size {} [e]'.format(label1))
        ax.set_ylabel('Jump Size {} [e]'.format(label2))
        ax.scatter(jumps1, jumps2, c=corrJumps, cmap='rainbow', marker='.')
        plt.draw()
        plt.pause(0.05)
        
    def plot_time_steps(self, datasets=None):
        if datasets is None:
            datasets = self.file_tags
        dt = np.array([])
        for ftag in datasets:
            dt = np.append(dt, self.time[ftag][1:] - self.time[ftag][:-1])
            
        print('Average Measurement Interval T = {:.2f}'.format(np.mean(dt)))
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Time Steps')
        ax.set_xlabel('N')
        ax.set_ylabel('Step Size [s]')
        ax.plot(dt)
        plt.draw()
        plt.pause(0.05)
        
        
