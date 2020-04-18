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
CO = ChargeOffset()
CO.add_dataset('fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvh0508hpk_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvn0611rgv_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvo0228okh_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvo0526mwl_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvs0514nqa_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvv0302hgc_correlation.hdf5')
CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                '\Parameter\cvw0030mes_correlation.hdf5')

files = []
path = ('Z:\\mcdermott-group\\data\\fluxNoise\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\03-20-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvh0508hpk', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\03-26-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvn0611rgv', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
# path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\03-26-20\\Charge_resetting\MATLABData\\')
# files += CO.get_files_triggered_on_bad_fit('cvo0228okh', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\03-27-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvo0526mwl', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\03-31-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvs0514nqa', 'Q{}'.format(Q_A), path.format('Q'+str(Q_A)))
path1 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\04-02-20\\Charge_resetting\MATLABData\\')
path2 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\04-03-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvv0302hgc', 'Q{}'.format(Q_A), [path1.format('Q'+str(Q_A)), path2.format('Q'+str(Q_A))] )
path1 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\04-03-20\\Charge_resetting\MATLABData\\')
path2 = ('Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}\\General\\04-04-20\\Charge_resetting\MATLABData\\')
files += CO.get_files_triggered_on_bad_fit('cvw0030mes', 'Q{}'.format(Q_A), [path1.format('Q'+str(Q_A)), path2.format('Q'+str(Q_A))] )
    
acfH = np.zeros(2*50+1)
acfL = np.zeros(2*50+1)
acf  = np.zeros(2*50+1)
acf0 = np.zeros(2*50+1)
avg_before = []
avg_after = []
avg_trial_before = []
avg_trial_after = []
avg_prev_file = []
add_P1_trace = np.zeros(2*10000+1)
add_flip_trace = np.zeros(2*10000+1)
n_trace = np.zeros(2*10000+1)

# whitelist = [2,7,16,28,44,46,51,52,61,67,80,106,114,124,131,145]
# whitelist = [2,7,16,44,46,51,52,61,106,114,124] # selected from average graph LOW
whitelist = [28,67,80] # selected from average graph HIGH
# whitelist = [2,7,16,28,44,51,52,61,106,114,124] # selected in zoom with robert by looking at events
blacklist = {74:range(10), 90:range(10), 92:range(10), 145:[3]}
data_files = {}

# for k,f in enumerate(files):
for k,f in [(i,files[i]) for i in whitelist]:
    path, num = noiselib.path_to_num(f)
    DF = TwoMeasDataFile(path.format(num+Q_A-1), 1, 2)
    DF.set_trigger_params(100000., 100, 2)
    # DF.apply_infidelity_correction(7)
    trigs = DF.get_triggers()
    trigs = [(t,r) for t,r in trigs if k not in blacklist.keys() 
                                    or t not in blacklist[k]]
    data_files[f] = DF
    print k, trigs
    for trial, rep in trigs:
        # DF.plot(trial)
        P1 = DF.get_P1_around_trig((trial,rep))
        nreps = P1.size
        n_trace[10000-rep:10000+nreps-rep] += 1
        add_P1_trace[10000-rep:10000+nreps-rep] += P1
        add_flip_trace[10000-rep:10000+nreps-rep-1] += np.abs(np.diff(P1))
        avg_before += [np.mean(P1[:rep])]
        avg_after += [np.mean(P1[rep:])]
        if trial > 0:
            P1 = DF.get_P1_around_trig((trial-1,rep))
            avg_trial_before += [np.mean(P1)]
        else:
            avg_trial_before += [np.nan]
        if trial < 9:
            P1 = DF.get_P1_around_trig((trial+1,rep))
            avg_trial_after += [np.mean(P1)]
        else:
            avg_trial_after += [np.nan]
        # avg_prev_file += [np.mean(P10)]
        # avg += [np.mean(P1)]
        
        i_0, i_end = max(0, rep-00), min(len(P1)-1, rep+100)
        acf  += noiselib.crosscorrelate(P1[i_0:i_end],  P1[i_0:i_end],  50)[0]
        # acf0 += noiselib.crosscorrelate(P10[i_0:i_end], P10[i_0:i_end], 50)[0]
        
        if np.mean(P1) > 0.06:
            acfH += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]
        else:
            acfL += noiselib.crosscorrelate(P1[i_0:i_end], P1[i_0:i_end], 50)[0]
    
# ax = plt.figure().add_subplot(111)
# ax.plot(np.arange(-50,51), acfH/np.sum(np.array(avg_after)>0.06), label='Triggered on Q{} charge jump high'.format(Q_A))
# ax.plot(np.arange(-50,51), acfL/np.sum(np.array(avg_after)<=0.06), label='Triggered on Q{} charge jump low'.format(Q_A))
# ax.plot(np.arange(-50,51), acf/len(avg_after), label='Triggered on Q{} charge jump all'.format(Q_A))
# ax.plot(np.arange(-50,51), acf0/len(avg_after), label='Triggered one file previous')
# ax.set_xlabel('Lag [100us]')
# ax.set_ylabel('Autocorrelation amplitude for Q{}'.format(Q_B))
# ax.legend()

ax = plt.figure().add_subplot(111)
ax.plot(add_P1_trace)
ax.plot(add_P1_trace/n_trace)
ax.plot(movingmean(add_P1_trace,100)/n_trace)
ax.plot(movingmean(add_P1_trace,1000)/n_trace)
ax.plot(n_trace)

ax = plt.figure().add_subplot(111)
ax.plot(avg_before, label='before')
ax.plot(avg_after, label='after')
ax.plot(avg_trial_before, label='previous')
ax.plot(avg_trial_after, label='next')
# ax.plot(avg_prev_file, label='prev file')
ax.legend()
ax.set_xticks(range(len(avg_before))) 
ax.set_xticklabels(whitelist)

plt.draw()
plt.pause(0.05)