import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from noiselib import movingmean

class TwoMeasDataFile(object):
    
    def __init__(self, path, charge_var='Single_Shot_Occupations_SB1', 
                             meas_var='Single_Shot_Occupations_SB2', meas_fn=None):
        self.path = path
        data = noiselib.loadmat( path )
        # self.o_charge = np.array(data[charge_var], dtype=np.bool)
        if isinstance(charge_var, str):
            self.o_charge = np.array(data[charge_var], dtype=np.bool)
        else:
            self.os_charge = [np.array(data[var], dtype=np.bool) for var in meas_var]
            self.o_charge = np.logical_xor(*self.os_charge)
        if isinstance(meas_var, str):
            self.o_meas = np.array(data[meas_var], dtype=np.bool)
            self.meas_var = meas_var
        else:
            self.os_meas = [np.array(data[var], dtype=np.bool) for var in meas_var]
            self.meas_var = meas_var[0]
            if meas_fn is None:
                self.o_meas = self.os_meas[0]
            else:
                self.o_meas = meas_fn(*self.os_meas)
    
    def apply_infidelity_correction(self, n_bins=9, thresh=0.5):
        self.o_meas = noiselib.apply_infidelity_correction(self.o_meas, n_bins, thresh)

    def apply_infidelity_correction_HMM(self, fidelity=[0.95, 0.75], trig=None, 
                                              fidelity2=[0.95, 0.75]):
        if trig is None:
            self.o_meas = noiselib.apply_infidelity_correction_HMM(self.o_meas, 
                                                                fidelity=fidelity)
        else:
            o1 = self.o_meas[trig[0]][:trig[1]]
            o2 = self.o_meas[trig[0]][trig[1]:]
            o1 = noiselib.apply_infidelity_correction_HMM( o1, fidelity=fidelity )
            o2 = noiselib.apply_infidelity_correction_HMM( o2, fidelity=fidelity2 )
            self.o_meas[trig[0]] = np.append(o1, o2)
    
    def add_base_ramsey_trace(self, data):
        self.base_trace = data

    def set_trigger_params(self, amp, span, thresh):
        self.trig = { 'amp': amp,
                      'span': span,
                      'thresh': thresh }
    
    def _trigger_trace(self, trial):
        o = self.o_charge[trial]
        smoothed_trial = movingmean(o, self.trig['span'])
        smoothed_gradient = movingmean(np.gradient(smoothed_trial), self.trig['span'] )
        return self.trig['amp']*smoothed_gradient**2
    
    def get_triggers(self, trials=None):
        trigs = []
        for _trial in range(len(self.o_charge)):
            if trials is not None and _trial not in trials:
                continue
            trig_trace = self._trigger_trace(_trial)
            rep = np.argmax(trig_trace)
            if trig_trace[rep] > self.trig['thresh']:
                trigs += [(_trial, rep)]
        return trigs
    
    def plot(self, trial, plot='P1', smoothing=20):
        path, num = noiselib.path_to_num(self.path)
        o = self.o_meas[trial]
        fig,axs = plt.subplots(4, 1, sharex=True)
        fig.suptitle('{}.{}'.format(num,trial))
        
        ax = axs[0]
        ax.get_shared_x_axes().remove(ax)
        if hasattr(self, 'base_trace'):
            ax.plot(self.base_trace, 'k:', label='')
        ax.plot(np.mean(self.o_charge, axis=1), label='')
        ax.plot([trial],[np.mean(self.o_charge, axis=1)[trial]], 'r.', markersize=18)
        ax.get_xaxis().set_ticks([])
        
        ax = axs[1]
        ax.plot(self._trigger_trace(trial), label='Trigger Trace')
        ax.plot(movingmean(self.o_charge[trial], 100), linewidth=0.5, 
                label='Smoothed P1 - Charge Qubit')
        ax.get_xaxis().set_ticks([])
        
        ax = axs[2]
        ax.plot(o, linewidth=0.2, label='P1')
        ax.get_xaxis().set_ticks([])
        
        ax = axs[3]
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
        return axs
        
    def plot_adj_file(self, path, x):
        if x == 'end':
            x = self.o_meas[0].size - 1
        data = noiselib.loadmat( path )
        o_meas = np.array(data['Single_Shot_Occupations_SB{}'.format(self.meas_var)], dtype=np.bool)
        o_meas = noiselib.apply_infidelity_correction(o_meas, 9, 0.5)
        means = np.mean(np.abs(np.diff(o_meas,axis=1)), axis=1)
        # self.ax.scatter( x, np.mean(means) )
        self.ax.errorbar( x, np.mean(means), np.std(means)/np.sqrt(len(means)), marker='.', mfc='green', mec='green', ecolor='green' )
    
    def get_P1_around_trig(self, trial_rep, reps_before=1e6, reps_after=1e6, meas_index=None):
        trial, rep = trial_rep
        if not meas_index:
            o = self.o_meas[trial]
        else:
            o = self.os_meas[meas_index][trial]
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