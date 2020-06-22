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
charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
               '\Parameter\{}_correlation.hdf5')
CO = ChargeOffset()
# # after T1 measurements, starting 6/13
# CO.add_dataset(charge_path.format('cyo0514vli'))
# CO.add_dataset(charge_path.format('cyp0520ynw'))
CO.add_dataset(charge_path.format('cyw0458jbl')).limit_dataset('Q1',start=77,end=122) \
                                                .limit_dataset('Q1',start=130,end=153) \
                                                .limit_dataset('Q1',start=205,end=243) \
                                                .limit_dataset('Q1',start=287,end=308) \
                                                .limit_dataset('Q1',start=378) \

files = []
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\Charge_resetting\MATLABData\\'
# # after T1 measurements, starting 6/13
# files += CO.get_files_triggered_on_bad_fit('cyo0514vli', 'Q{}'.format(Q_A), path.format('06-13-20'))
# files += CO.get_files_triggered_on_bad_fit('cyp0520ynw', 'Q{}'.format(Q_A), path.format('06-14-20'))
files += CO.get_files_triggered_on_bad_fit('cyw0458jbl', 'Q{}'.format(Q_A), [path.format('06-20-20'),path.format('06-21-20')])

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

whitelist = None
# whitelist = {}
blacklist = {}
data_files = {}

print len(files)
for k,f in enumerate(files[0:50]):
    path, num = noiselib.path_to_num(f)
    date = f.split('\\')[-4]
    if (whitelist is not None 
        and (date not in whitelist or num not in whitelist[date])):
        continue
    DF = TwoMeasDataFile(path.format(num+Q_A), charge_var='Single_Shot_Occupations_RO2_SB1', 
                                                 meas_var=['Single_Shot_Occupations_RO2_SB2',
                                                           'Single_Shot_Occupations_RO1_SB2'])
    DF.set_trigger_params(1000000., 300, 0.5)
    DF.apply_infidelity_correction_HMM(fidelity=[0.95, 0.75])
    # DF.apply_infidelity_correction(9)
    trigs = DF.get_triggers()
    trigs = [(t,r) for t,r in trigs if k not in blacklist.keys() 
                                    or t not in blacklist[k]]
    if len(trigs) > 0:
        data_files[f] = DF
    print k, trigs
    for trial, rep in trigs:
        DF.plot(trial, plot='Flips', smoothing=500)
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

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(np.arange(-50,51), acfH/np.sum(np.array(avg)>0.06), label='Triggered on Q{} charge jump high'.format(Q_A))
# ax.plot(np.arange(-50,51), acfL/np.sum(np.array(avg)<=0.06), label='Triggered on Q{} charge jump low'.format(Q_A))
# ax.plot(np.arange(-50,51), acf/len(avg), label='Triggered on Q{} charge jump all'.format(Q_A))
# ax.plot(np.arange(-50,51), acf0/len(avg), label='Triggered one file previous')
# ax.set_xlabel('Lag [100us]')
# ax.set_ylabel('Autocorrelation amplitude for Q{}'.format(Q_B))
# ax.legend()
# plt.draw()
# plt.pause(0.05)

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
ax.set_xlabel('Rep (trigger centered on 10k)')
ax.set_ylabel('Smoothed Flip Rate')
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
