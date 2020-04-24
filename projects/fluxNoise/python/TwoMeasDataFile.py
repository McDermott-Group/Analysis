import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean

class TwoMeasDataFile(object):
    
    def __init__(self, path, charge_SB, meas_SB):
        self.path = path
        self.meas_SB = meas_SB
        data = noiselib.loadmat( path )
        self.o_charge = np.array(data['Single_Shot_Occupations_SB{}'.format(charge_SB)], dtype=np.bool)
        self.o_meas = np.array(data['Single_Shot_Occupations_SB{}'.format(meas_SB)], dtype=np.bool)
    
    def apply_infidelity_correction(self, n_bins=9, thresh=0.5):
        self.o_meas = noiselib.apply_infidelity_correction(self.o_meas, n_bins, thresh)
    
    def set_trigger_params(self, amp, span, thresh):
        self.trig = { 'amp': amp,
                      'span': span,
                      'thresh': thresh }
    
    def _trigger_trace(self, trial):
        o = self.o_charge[trial]
        smoothed_trial = movingmean(o, self.trig['span'])
        smoothed_gradient = movingmean(np.gradient(smoothed_trial), self.trig['span'] )
        return self.trig['amp']*smoothed_gradient**2
    
    def get_triggers(self):
        trigs = []
        for trial in range(len(self.o_charge)):
            trig_trace = self._trigger_trace(trial)
            rep = np.argmax(trig_trace)
            if trig_trace[rep] > self.trig['thresh']:
                trigs += [(trial, rep)]
        return trigs
    
    def plot(self, trial, plot='P1', smoothing=20):
        path, num = noiselib.path_to_num(self.path)
        o = self.o_meas[trial]
        fig = plt.figure()
        
        ax = fig.add_subplot(311)
        ax.set_title(num)
        
        ax.plot(self._trigger_trace(trial), label='Trigger Trace')
        ax.plot(movingmean(self.o_charge[trial], 100), linewidth=0.5, label='Smoothed P1 - Charge Qubit')
        ax.get_xaxis().set_ticks([])
        
        ax = fig.add_subplot(312)
        ax.plot(o, linewidth=0.2, label='P1')
        ax.get_xaxis().set_ticks([])
        
        ax = fig.add_subplot(313)
        end = self.o_meas[0].size - 1
        if plot == 'P1':
            ax.plot(movingmean(o, smoothing), label='Smoothed P1')
            if trial > 0:
                adjacent = np.mean(self.o_meas[trial-1])
                ax.plot([0, 1000], [adjacent]*2, 'ro-', markevery=[0], markersize=5)
                ax.plot([end-1000, end], [adjacent]*2, 'r--')
            if trial < len(self.o_meas)-1:
                adjacent = np.mean(self.o_meas[trial+1])
                ax.plot([0, 1000], [adjacent]*2, 'r--')
                ax.plot([end-1000, end], [adjacent]*2, 'ro-', markevery=[1], markersize=5)
        if plot == 'Flips':
            ax.plot(movingmean(np.abs(np.diff(o)), smoothing), label='Smoothed Flips')
            if trial > 0:
                adjacent = np.mean(np.abs(np.diff(self.o_meas[trial-1])))
                # ax.plot([0, 1000], [adjacent]*2, 'ro-', markevery=[0], markersize=5)
                ax.plot([end-1000, end], [adjacent]*2, 'r--')
                means_before = np.mean(np.abs(np.diff(self.o_meas[:trial],axis=1)), axis=1)
                ax.plot(100*np.arange(len(means_before)), means_before, 'r.-')
            if trial < len(self.o_meas)-1:
                adjacent = np.mean(np.abs(np.diff(self.o_meas[trial+1])))
                ax.plot([0, 1000], [adjacent]*2, 'r--')
                # ax.plot([end-1000, end], [adjacent]*2, 'ro-', markevery=[1], markersize=5)
                means_after = np.mean(np.abs(np.diff(self.o_meas[trial+1:],axis=1)), axis=1)
                ax.plot(end-100*np.arange(len(means_after)-1,-1,-1), means_after, 'r.-')
            
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        self.ax = ax
        return ax
        
    def plot_adj_file(self, path, x):
        if x == 'end':
            x = self.o_meas[0].size - 1
        data = noiselib.loadmat( path )
        o_meas = np.array(data['Single_Shot_Occupations_SB{}'.format(self.meas_SB)], dtype=np.bool)
        o_meas = noiselib.apply_infidelity_correction(o_meas, 9, 0.5)
        means = np.mean(np.abs(np.diff(o_meas,axis=1)), axis=1)
        # self.ax.scatter( x, np.mean(means) )
        self.ax.errorbar( x, np.mean(means), np.std(means)/np.sqrt(len(means)), marker='.', mfc='green', mec='green', ecolor='green' )
    
    def get_P1_around_trig(self, trial_rep, reps_before=1e6, reps_after=1e6):
        trial, rep = trial_rep
        o = self.o_meas[trial]
        i_0 = max(0, rep - reps_before)
        i_max = min(o.size, rep + reps_after)
        return o[i_0:i_max]
        
    def get_averages(self):
        pass
    
    # def find_jump_in_ramsey(o1, o2, plot=False):
        # badIndicies = np.argwhere(np.abs(o1 - o2) > 0.1)
        # badIndex = np.maximum([0],badIndicies[0]-1) if len(badIndicies) > 0 else [0]
        # if plot:
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # ax.plot(o1)
            # ax.plot(o2)
            # ax.plot([badIndex],[o1[badIndex]], 'ro')
            # plt.draw()
            # plt.pause(0.05)
        # return badIndex