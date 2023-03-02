import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from noiselib import movingmean
from QPTunneling import *
import ChargeOffset
importlib.reload(ChargeOffset)
from ChargeOffset import *
import TwoMeasDataFile
importlib.reload(TwoMeasDataFile)
from TwoMeasDataFile import *
from dataChest import dataChest
from random import randrange
import datasets
importlib.reload(datasets)
import datasets as ds
import shutil


# Q_A, Q_B = 1,2
# CO = ChargeOffset()
# files_QP = []
# files_XIY_upper = []
# files_XIY_lower = []
# charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
               # '\Parameter\{}_correlation.hdf5')
# path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\{}\MATLABData\\'

# # for d in [ds.q1q2_0701_charge_QP, ds.q1q2_0702_charge_QP]: # pre-calibration only
# for d in [ds.q1q2_0703_charge_QP, ds.q1q2_0704_charge_QP, ds.q1q2_0706_charge_QP]: # pre,post cal, 500us
# # for d in [ds.q1q2_0707_charge_QP, ds.q1q2_0708_charge_QP]: # 30us
# # for d in [ds.q1q2_0708_charge_QP]: #new
          
    # CO.add_dataset(charge_path.format(d['path_charge']))
    # if 'exclude_data' in d:
        # for start,end in d['exclude_data']:
            # CO.limit_dataset('Q1',start=start,end=end)
        
    # files_QP += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                # [path.format(date,'Charge_resetting') for date in d['date']], d['thresh_charge'])
    # files_XIY_upper += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                # [path.format(date,'Fidelity_UpperBand') for date in d['date']], d['thresh_charge'])
    # files_XIY_lower += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                # [path.format(date,'Fidelity_LowerBand') for date in d['date']], d['thresh_charge'])
    

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

def change_file_num(f, delta):
    path, num = noiselib.path_to_num(f)
    return path.format(num + delta)

whitelist = None
# whitelist = {'07-02-20': [1.4, 27.1, 113.2, 113.7, 197.7, 207.1, 263.6, 271.2, 
                          # 331.0, 593.6],  # 161.8
             # '07-03-20': [94.4, 130.9, 288.5, 358.6, 400.4, 406.6, 408.8, 446.1, 
                          # 522.4, 672.0, 710.6, 720.1, 757.7, 869.8, 879.6, 883.5, 
                          # 899.6, 901.5, 981.0, 1017.4, 1025.8],
             # '07-04-20': [4.4, 24.5, 42.1, 50.0, 116.0, 138.7, 168.7, 310.0, 
                          # 362.5, 370.0, 388.3, 414.7, 430.5, 452.9, 460.5, 
                          # 480.4, 552.0, 610.7, 718.1, 822.1, 876.3, 886.3, 906.2],
             # '07-05-20': [134.1, 222.5, 250.4, 252.6, 268.8, 322.0, 328.2, 
                          # 420.5, 446.8, 906.2],
             # '07-06-20': [55.0, 179.8, 305.6, 403.9, 431.5,  957.4, 981.0],
             # '07-07-20': [],
             # '07-08-20': [499.6, 515.8, 933.1],
             # '07-09-20': [379.6, 623.7] 
             # }
blacklist = {}
data_files = {}

print(len(files_QP))
for k,f in enumerate(files_QP[10:70]):
    path, num = noiselib.path_to_num(f)
    num = num + Q_A
    date = f.split('\\')[-4]
    if (whitelist is not None 
        and (date not in whitelist or num not in np.floor(whitelist[date]))):
        continue
    DF = TwoMeasDataFile(path.format(num), charge_var='Single_Shot_Occupations_RO2_SB1', 
                                             meas_var=['Single_Shot_Occupations_RO2_SB2',
                                                       'Single_Shot_Occupations_RO1_SB2'])
    DF.set_trigger_params(1e7, 500, 1)
    trials = None
    if whitelist is not None and np.array(whitelist[date]).dtype == float:
        i = np.argwhere( np.floor(whitelist[date]) == num )
        v = np.array(whitelist[date])[i]
        trials = np.round((10*(v%1))).reshape(-1).astype(int)
    trigs = DF.get_triggers(trials)
    trigs = [(t,r) for t,r in trigs if date not in blacklist
                                    or num not in blacklist[date]
                                    or t not in blacklist[date][num]]
    
    
    if len(trigs) > 0:
        # data_files[f] = DF
        
        # s = 'Single_Shot_Occupation_RO2_SB2'
        # data0 = noiselib.loadmat( change_file_num(files_XIY_upper[k],-1) )
        # data = noiselib.loadmat( files_XIY_upper[k] )
        # XIY_upper = (data0[s], data[s])
        # data0 = noiselib.loadmat( change_file_num(files_XIY_lower[k],-1) )
        # data = noiselib.loadmat( files_XIY_lower[k] )
        # XIY_lower = (data0[s], data[s])
        # # print XIY_upper, XIY_lower
        
        s1 = 'Single_Shot_Occupations_RO1_SB2'
        s2 = 'Single_Shot_Occupations_RO2_SB2'
        data0 = noiselib.loadmat( change_file_num(files_XIY_upper[k],-1) )
        data = noiselib.loadmat( files_XIY_upper[k] )
        XIY_upper = (1.*np.sum(data0[s1]!=data0[s2])/data0[s2].size, 
                     1.*np.sum(data[s1]!=data[s2])/data[s2].size)
        data0 = noiselib.loadmat( change_file_num(files_XIY_lower[k],-1) )
        data = noiselib.loadmat( files_XIY_lower[k] )
        XIY_lower = (1.*np.sum(data0[s1]!=data0[s2])/data0[s2].size, 
                     1.*np.sum(data[s1]!=data[s2])/data[s2].size)
        # print XIY_upper, XIY_lower
        # print
        
    print(k, date, num, trigs)
    for trial, rep in trigs:
        # copy_path = 'C:\\Users\\lab_user\\Desktop\\for_livermore'
        # shutil.copy2( path.format(num), copy_path )
        # shutil.copy2( change_file_num(files_XIY_upper[k],-1), copy_path )
        # shutil.copy2( files_XIY_upper[k], copy_path )
        # shutil.copy2( change_file_num(files_XIY_lower[k],-1), copy_path )
        # shutil.copy2( files_XIY_lower[k], copy_path )
        # print path.format(num)
        # print '\t', change_file_num(files_XIY_upper[k],-1)
        # print '\t', files_XIY_upper[k]
        # print '\t', change_file_num(files_XIY_lower[k],-1)
        # print '\t', files_XIY_lower[k]
        # print
        
        DF.apply_infidelity_correction_HMM(fidelity=[XIY_upper[0], XIY_upper[0]],
                                           trig=(trial,rep),
                                           fidelity2=[XIY_upper[1], XIY_upper[1]])
        # DF.apply_infidelity_correction_HMM(fidelity=[0.95, 0.85])
        axs = DF.plot(trial, plot='Flips', smoothing=500)
        # DF.plot_adj_file(path.format(num-3), 0)
        # DF.plot_adj_file(path.format(num+1), 'end')
        P1 = DF.get_P1_around_trig((trial,rep))
        nreps = P1.size
        axs[1].text(0, 0.3, '{:.2f}\n{:.2f}'.format(XIY_lower[0], XIY_upper[0]))
        axs[1].text(nreps, 0.3, '{:.2f}\n{:.2f}'.format(XIY_lower[1], XIY_upper[1]), horizontalalignment='right')
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
