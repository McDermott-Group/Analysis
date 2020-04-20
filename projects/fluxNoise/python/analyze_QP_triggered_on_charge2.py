import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
import TwoMeasDataFile
reload(TwoMeasDataFile)
from TwoMeasDataFile import *
from dataChest import dataChest
from random import randrange

Q_A, Q_B = 1,2
offset_path = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter'
CO = ChargeOffset()
CO.add_dataset(offset_path + '\cvi1713bxq_correlation.hdf5')
CO.add_dataset(offset_path + '\cvo2209zfd_correlation.hdf5')
CO.add_dataset(offset_path + '\cvp0501ixo_correlation.hdf5')
CO.add_dataset(offset_path + '\cvt0001lvs_correlation.hdf5')

files = []
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\Charge_resetting_QP\MATLABData\\'
# files += CO.get_files_triggered_on_bad_fit('cvi1713bxq', 'Q{}'.format(Q_A), path.format('03-21-20'))
files += CO.get_files_triggered_on_bad_fit('cvo2209zfd', 'Q{}'.format(Q_A), path.format('03-27-20'))
files += CO.get_files_triggered_on_bad_fit('cvp0501ixo', 'Q{}'.format(Q_A), path.format('03-28-20'))
files += CO.get_files_triggered_on_bad_fit('cvt0001lvs', 'Q{}'.format(Q_A), [path.format('03-31-20'), 
                                                                             path.format('04-01-20')])
    
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
add_trial_means = np.zeros(2*10+1)
n_trials = np.zeros(2*10+1)

whitelist = [2,3,4,12,38,40,45,48,51,52,54,71,74,75,78,79,86]
blacklist = {40:[7], 34:[4], 27:[7,3], 26:[3,1], 24:[7], 20:[8], 19:[7,5,2,0], 18:[8], 17:[9,3], 13:[7,5], 11:[9,3,1], 9:[3], 8:[4], 6:[5,3], 3:[8,6,3], 0:[6,4]}
data_files = {}

# for k,f in enumerate(files):
for k,f in [(i,files[i]) for i in whitelist[:]]:
    path, num = noiselib.path_to_num(f)
    DF = TwoMeasDataFile(path.format(num+Q_A-1), 1, 2)
    DF.set_trigger_params(1000000., 1000, 0.05)
    DF.apply_infidelity_correction(9)
    trigs = DF.get_triggers()
    trigs = [(t,r) for t,r in trigs if k not in blacklist.keys() 
                                    or t not in blacklist[k]]
    if len(trigs) > 0:
        data_files[f] = DF
    print k, trigs
    for trial, rep in trigs:
        # DF.plot(trial, plot='Flips', smoothing=500)
        # DF.plot_adj_file(path.format(num+Q_A-3), 0)
        # DF.plot_adj_file(path.format(num+Q_A+1), 'end')
        P1 = DF.get_P1_around_trig((trial,rep))
        nreps = P1.size
        n_trace[10000-rep:10000+nreps-rep] += 1
        add_P1_trace[10000-rep:10000+nreps-rep] += P1
        add_flip_trace[10000-rep:10000+nreps-rep-1] += np.abs(np.diff(P1))
        avg_before += [np.mean(P1[:rep])]
        avg_after += [np.mean(P1[rep:])]
        add_trial_means[10-trial:20-trial] += np.mean(np.abs(np.diff(DF.o_meas, axis=1)), axis=1)
        n_trials[10-trial:20-trial] += 1
        # avg_prev_file += [np.mean(P10)]
        avg += [np.mean(P1)]
        
        i_0, i_end = max(0, rep-00), min(len(P1)-1, rep+100)
        acf  += noiselib.crosscorrelate(P1[i_0:i_end],  P1[i_0:i_end],  50)[0]
        # acf0 += noiselib.crosscorrelate(P10[i_0:i_end], P10[i_0:i_end], 50)[0]
        
        if np.mean(P1) > 0.06:
            acfH += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]
        else:
            acfL += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]

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

# ax = plt.figure().add_subplot(111)
# ax.plot(add_P1_trace)
# avg_P1 = np.zeros(2*10000+1)
# ax.plot(np.divide(add_P1_trace, n_trace, avg_P1, where=n_trace>0))
# ax.plot(n_trace)

ax = plt.figure().add_subplot(111)
ax.plot(add_flip_trace)
avg_flip = np.zeros(2*10000+1)
avg_flip_trace = np.divide(add_flip_trace, n_trace, avg_flip, where=n_trace>0)
ax.plot(avg_flip_trace)
ax.plot(movingmean(avg_flip_trace, 100))
ax.plot(movingmean(avg_flip_trace, 1000))
ax.plot(n_trace)
plt.draw()
plt.pause(0.05)

ax = plt.figure().add_subplot(111)
trial_avgs = np.zeros(2*10+1)
trial_avgs = np.divide(add_trial_means, n_trials, trial_avgs, where=n_trials>0)
ax.plot(trial_avgs)
ax.plot(n_trials)
ax.set_xlabel('Trial')
ax.set_ylabel('Mean Flip Rate')
plt.draw()
plt.pause(0.05)


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