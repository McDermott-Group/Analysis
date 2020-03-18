import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt

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
                offset[l] = np.append( offset[l], offset[l][-1] + new_charge )
            for l in set(offset.keys()) - set(self.labels[ftag]):
                offset[ftag] = np.append( offset[ftag], np.full(n_new, np.nan) )
        return time, offset
    
    def get_jump_sizes(self, datasets = None):
        time, offset = self.get_charge_offset(datasets)
        return {l: offset[l][1:] - offset[l][:-1] for l in offset.keys()}
        
    def plot_charge_offset(self, datasets = None):
        time, offset = self.get_charge_offset(datasets)
        jumps = self.get_jump_sizes(datasets)
        largeJumps = {l: np.append(abs(jumps[l]) > 0.1, True) | np.append(True, abs(jumps[l]) > 0.1)
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
        jumps = self.get_jump_sizes(datasets)
        for l in jumps.keys():
            n_tot = jumps[l].size
            n_big = np.sum(abs(jumps[l]) > 0.1)
            print('{}: {} / {} = {:.2f}% = {:.2f}*T between jumps'.format(
                    l, n_big, n_tot, 100.*n_big/n_tot, 1.*n_tot/n_big))
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Jump Sizes')
        ax.set_xlabel('N')
        ax.set_ylabel('Jump Size [e]')
        for l in jumps.keys():
            ax.plot(np.abs(jumps[l]), label=l)
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        
    def plot_charge_correlation(self, label1, label2, datasets=None):
        jumps = self.get_jump_sizes(datasets)
        jumps1 = jumps[label1]
        jumps2 = jumps[label2]
        corrJumps = (abs(jumps1) > 0.1) & (abs(jumps2) > 0.1)
        
        print( 'Correlated Jumps: {}-{}'.format(label1,label2) )
        for j in [jumps1,jumps2]:
            print( '    {} / {} = {:.2f}%'.format( 
                    np.sum(corrJumps), 
                    np.sum(abs(j) > 0.1),
                    100.*np.sum(corrJumps)/np.sum(abs(j) > 0.1) ) )
                
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
        
        
base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Q3Q4Corr\General\Parameter\\'
CO = ChargeOffset()
# CO.add_dataset(base_path + 'cvc0232imo_correlation.hdf5') # bad Q2
# CO.add_dataset(base_path + 'cvd0430bxy_correlation.hdf5') # bad Q4
# CO.add_dataset(base_path + 'cvd0802ltu_correlation.hdf5') # bad Q4
CO.add_dataset(base_path + 'cvd1745asc_correlation.hdf5')
CO.add_dataset(base_path + 'cve0104ppt_correlation.hdf5')
CO.add_dataset(base_path + 'cvf0542fei_correlation.hdf5')
CO.plot_charge_offset()
CO.plot_jump_sizes()
CO.plot_charge_correlation('Q1','Q2')
CO.plot_charge_correlation('Q3','Q4')
CO.plot_charge_correlation('Q1','Q3')
CO.plot_time_steps()
plt.show()

"""######### Interleaved Data ##########"""

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q3Q4Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csn0100hzy_correlation.hdf5')
# CO.add_dataset(base_path + 'csn1803vvg_correlation.hdf5')
# CO.add_dataset(base_path + 'cso0030jjg_correlation.hdf5')
# CO.add_dataset(base_path + 'csq0128zge_correlation.hdf5')
# # CO.add_dataset(base_path + 'csy2025oeg_correlation.hdf5')
# # CO.add_dataset(base_path + 'csz0145znf_correlation.hdf5')
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q3','Q4')
# CO.plot_time_steps()
# plt.show()

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csu0052mrt_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0034shs_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0824hxv_correlation.hdf5')
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q2')
# CO.plot_time_steps()
# plt.show()

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q3Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csw0353opr_correlation.hdf5')
# CO.add_dataset(base_path + 'csx0545cpf_correlation.hdf5')
# CO.add_dataset(base_path + 'csx1710hcz_correlation.hdf5')
# CO.add_dataset(base_path + 'csy0726yvn_correlation.hdf5')
# # base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q4\General\Parameter\\'
# # CO.add_dataset(base_path + 'crc2345hmc_parameters.hdf5') # might say Q3,Q4 instead?
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q3')
# CO.plot_time_steps()
# plt.show()