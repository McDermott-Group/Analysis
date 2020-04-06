import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
from dataChest import dataChest
from random import randrange

def find_triggers_in_file(o, threshold=2, o2=None, plot=False, title=''):
    trigs = []
    for i,trial in enumerate(o):
        trig_trace = 1000000.*movingmean(np.gradient(movingmean(trial, 1000)), 1000 )**2
        if trig_trace.max() > threshold: # and trig_trace.argmax() > 300 and trig_trace.argmax() < 9700:
            trigs += [(i, trig_trace.argmax())]
            if plot:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title(title)
                ax.plot(trig_trace, label='Trigger Trace')
                ax.plot(movingmean(trial, 100), label='Smoothed P1 - Charge Qubit')
                if o2 is not None:
                    P1_B = o2[i]
                    ax.plot(-3 + 0.5*o2[i], linewidth=0.2, label='P1 - QP Qubit')
                    ax.plot(-2 + movingmean(o2[i], 100), label='Smoothed P1 - QP Qubit')
                    ax.plot(-1 + 2*movingmean(np.abs(np.gradient(o2[i])), 100), label='Smoothed Gradient - QP Qubit')
                plt.legend()
                plt.draw()
                plt.pause(0.05)
    return trigs

def find_jump_in_ramsey(o1, o2, plot=False):
    badIndicies = np.argwhere(np.abs(o1 - o2) > 0.1)
    badIndex = np.maximum([0],badIndicies[0]-1) if len(badIndicies) > 0 else [0]
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(o1)
        ax.plot(o2)
        ax.plot([badIndex],[o1[badIndex]], 'ro')
        plt.draw()
        plt.pause(0.05)
    return badIndex


# class DataFile(object):
    
    # def __init__(self, path, charge_SB, meas_SB):
        # self.path = path
        # data = noiselib.loadmat( path )
        # self.o_charge = np.array(data['Single_Shot_Occupations_SB{}'.format(charge_SB)], dtype=np.bool)
        # self.o_meas = np.array(data['Single_Shot_Occupations_SB{}'.format(meas_SB)], dtype=np.bool)
    
    # def apply_infidelity_correction(self, n_bins=7, thresh=0.5):
        # for trial in self.o_meas:
            # trial[...] = movingmean(trial, n_bins) > thresh
    
    # def set_trigger_params(self, amp, span, thresh):
        # self.trig = { 'amp': amp,
                      # 'span': span,
                      # 'thresh': thresh }
    
    # def _trigger_trace(self, trial):
        # o = self.o_charge[trial]
        # smoothed_trial = movingmean(o, self.trig['span'])
        # smoothed_gradient = movingmean(np.gradient(smoothed_trial), self.trig['span'] )
        # return self.trig['amp']*smoothed_gradient**2
    
    # def get_triggers(self):
        # trigs = []
        # for trial in range(len(self.o_charge)):
            # trig_trace = self._trigger_trace(trial)
            # rep = np.argmax(trig_trace)
            # if trig_trace[rep] > self.trig['thresh']:
                # trigs += [(trial, rep)]
        # return trigs
    
    # def plot(self, trial):
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # path, num = noiselib.path_to_num(f)
        # ax.set_title(num)
        
        # o = self.o_meas[trial]
        # ax.plot(1 + self._trigger_trace(trial), label='Trigger Trace')
        # ax.plot(1 + movingmean(self.o_charge[trial], 100), linewidth=0.5, label='Smoothed P1 - Charge Qubit')
        # ax.plot(0.5*o, linewidth=0.2, label='P1')
        # ax.plot(-1 + movingmean(o, 100), label='Smoothed P1')
        # ax.plot(-2 + 2*movingmean(np.abs(np.gradient(o)), 100), label='Smoothed Gradient')
        # avg_adjacent = []
        # if trial > 0:
            # avg_adjacent += [(0,np.mean(self.o_meas[trial-1]))]
        # if trial < len(self.o_meas-1):
            # avg_adjacent += [(self.o_meas[0].size,np.mean(self.o_meas[trial+1]))]
        # ax.scatter(*zip(*avg_adjacent), 'o')
            
        # ax.legend()
        # plt.draw()
        # plt.pause(0.05)
        # return ax
    
    # def get_P1_around_trig(self, trial_rep, reps_before=None, reps_after=None):
        # pass
        
    # def get_averages(self):
        # pass

Q_A, Q_B = 1,2
CO = ChargeOffset()
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvi1713bxq_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvo2209zfd_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvp0501ixo_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvt0001lvs_correlation.hdf5')

files = []
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\03-21-20\\Charge_resetting_QP\MATLABData\\')
# files += CO.get_files_triggered_on_bad_fit('cvi1713bxq', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\03-27-20\\Charge_resetting_QP\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvo2209zfd', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\03-28-20\\Charge_resetting_QP\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvp0501ixo', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path1 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\03-31-20\\Charge_resetting_QP\MATLABData\\')
path2 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\04-01-20\\Charge_resetting_QP\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvt0001lvs', 'Q{}'.format(Q_A), [path1.format('Q'+str(Q_A)), path2.format('Q'+str(Q_A))])
    
acfH = np.zeros(2*50+1)
acfL = np.zeros(2*50+1)
acf  = np.zeros(2*50+1)
acf0 = np.zeros(2*50+1)
avg = []
avg_before = []
avg_after = []
avg_prev_file = []
add_P1_trace = np.zeros(2*10000+1)
add_flip_trace = np.zeros(2*10000+1)
n_trace = np.zeros(2*10000+1)

whitelist = [2,3,4,12,38,40,45,48,51,52,54,71,74,75,78,79,86]
blacklist = {40:[7], 34:[4], 27:[7,3], 26:[3,1], 24:[7], 20:[8], 19:[7,5,2,0], 18:[8], 17:[9,3], 13:[7,5], 11:[9,3,1], 9:[3], 8:[4], 6:[5,3], 3:[8,6,3], 0:[6,4]}

# for k,f in enumerate(files):
for k,f in [(i,files[i]) for i in whitelist[14:]]:
    path, num = noiselib.path_to_num(f)
    data = noiselib.loadmat(path.format(num+Q_A-1))
    SBA = np.array(data['Single_Shot_Occupation_SB{}'.format(Q_A)])
    SBAs = np.array(data['Single_Shot_Occupations_SB{}'.format(Q_A)])
    SBBs = np.array(data['Single_Shot_Occupations_SB{}'.format(Q_B)])
    for trial in SBBs:
        trial[...] = movingmean(trial, 7) > 0.5
    # data0 = noiselib.loadmat(path.format(num+Q_A-3))
    # SBBs0 = np.array(data0['Single_Shot_Occupations_SB{}'.format(Q_B)])
    # find_jump_in_ramsey(SBA, SBA0, plot=True)
    trigs = find_triggers_in_file(SBAs, .05, SBBs, plot=False, title=str(k))
    trigs = [(t,r) for t,r in trigs if k not in blacklist.keys() 
                                    or t not in blacklist[k]]
    print k, trigs,
    for trial, rep in trigs:
        P1 = SBBs[trial]
        P1 = movingmean(P1, 7) > 0.5
        nreps = P1.size
        # P10 = SBBs0[trial]
        n_trace[10000-rep:10000+nreps-rep] += 1
        add_P1_trace[10000-rep:10000+nreps-rep] += P1
        add_flip_trace[10000-rep:10000+nreps-rep] += np.abs(np.gradient(P1))
        avg_before += [np.mean(P1[:rep])]
        avg_after += [np.mean(P1[rep:])]
        # avg_prev_file += [np.mean(P10)]
        avg += [np.mean(P1)]
        
        i_0, i_end = max(0, rep-00), min(len(P1)-1, rep+100)
        acf  += noiselib.crosscorrelate(P1[i_0:i_end],  P1[i_0:i_end],  50)[0]
        # acf0 += noiselib.crosscorrelate(P10[i_0:i_end], P10[i_0:i_end], 50)[0]
        
        if np.mean(P1) > 0.06:
            acfH += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]
        else:
            acfL += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]
        # plt.figure()
        # plt.plot(P1, linewidth=0.2)
        # plt.plot(1 + 2*movingmean(P1, 20) )
        # plt.plot(2 + 2*movingmean(np.abs(np.gradient(P1)), 20) )
        # plt.draw()
        # plt.pause(0.05)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.arange(-50,51), acfH/np.sum(np.array(avg)>0.06), label='Triggered on Q{} charge jump high'.format(Q_A))
ax.plot(np.arange(-50,51), acfL/np.sum(np.array(avg)<=0.06), label='Triggered on Q{} charge jump low'.format(Q_A))
ax.plot(np.arange(-50,51), acf/len(avg), label='Triggered on Q{} charge jump all'.format(Q_A))
ax.plot(np.arange(-50,51), acf0/len(avg), label='Triggered one file previous')
ax.set_xlabel('Lag [100us]')
ax.set_ylabel('Autocorrelation amplitude for Q{}'.format(Q_B))
ax.legend()
plt.draw()
plt.pause(0.05)

ax = plt.figure().add_subplot(111)
ax.plot(add_P1_trace)
avg_P1 = np.zeros(2*10000+1)
ax.plot(np.divide(add_P1_trace, n_trace, avg_P1, where=n_trace>0))
ax.plot(n_trace)

ax = plt.figure().add_subplot(111)
ax.plot(add_flip_trace)
avg_flip = np.zeros(2*10000+1)
avg_flip_trace = np.divide(add_flip_trace, n_trace, avg_flip, where=n_trace>0)
ax.plot(avg_flip_trace)
ax.plot(movingmean(avg_flip_trace, 10))
ax.plot(n_trace)


# def bin_dwell_times(o):
     # last_val = o[0]
     # count = 0
     # x0 = []
     # x1 = []
     # for v in o:
         # if v == last_val:
             # count += 1
         # else:
             # if last_val == 0:
                 # x0 += [count]
             # else:
                 # x1 += [count]
             # count = 1
         # last_val = v
     # return x0,x1