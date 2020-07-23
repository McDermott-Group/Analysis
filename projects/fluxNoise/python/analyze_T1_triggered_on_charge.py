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
import datasets
reload(datasets)
import datasets as ds

Q_A, Q_B = 1,2
CO = ChargeOffset()
files = []
charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
               '\Parameter\{}_correlation.hdf5')
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\{}\MATLABData\\'

# for d in [ds.q1q2_0524_charge_T1, ds.q1q2_0524_charge_T1_2, ds.q1q2_0525_charge_T1, 
          # ds.q1q2_0526_charge_T1, ds.q1q2_0526_charge_T1_2, ds.q1q2_0527_charge_T1, 
          # ds.q1q2_0528_charge_T1, ds.q1q2_0528_charge_T1_2]: # 200us cycle, 3us wait
# for d in [ds.q1q2_0528_charge_T1_3, ds.q1q2_0529_charge_T1, ds.q1q2_0530_charge_T1, 
          # ds.q1q2_0603_charge_T1]: # 500us cycle, 3us wait
# for d in [ds.q1q2_0603_charge_T1_2, ds.q1q2_0604_charge_T1, ds.q1q2_0605_charge_T1,
          # ds.q1q2_0607_charge_T1, ds.q1q2_0608_charge_T1]: # 500us cycle, 9us wait
for d in [ds.q1q2_0711_charge_T1, ds.q1q2_0714_charge_T1, ds.q1q2_0715_charge_T1,
            ds.q1q2_0719_charge_T1, ds.q1q2_0720_charge_T1, ds.q1q2_0722_charge_T1]: # 30us duty cycle, 10us wait
# for d in [ds.q1q2_0722_charge_T1]: #new
          
    CO.add_dataset(charge_path.format(d['path_charge']))
    if 'exclude_data' in d:
        for start,end in d['exclude_data']:
            CO.limit_dataset('Q1',start=start,end=end)
        
    files += CO.get_files_triggered_on_bad_fit( d['path_charge'], 'Q{}'.format(Q_A), 
                [path.format(date,'Charge_resetting') for date in d['date']], 
                d['thresh_charge'])

    
N = 10000
P10 = np.zeros(2*N+1)
P01 = np.zeros(2*N+1)
P00 = np.zeros(2*N+1)
P11 = np.zeros(2*N+1)
n_trace = np.zeros(2*N+1)
n0 = np.zeros(2*N+1)
n1 = np.zeros(2*N+1)


whitelist = None
##### all events ######
# whitelist = {'05-24-20':[193,289,422,478,1347,1510],
             # '05-26-20':[560],
             # '05-27-20':[411],
             # '05-28-20':[1445,1705,1929],#1929
             # '05-29-20':[13,495,544,576,1514,1569,1867], #325
             # '05-30-20':[34,535,787],
             # '05-31-20':[92,317,1033,1079,1151,1250,1356,1539,1890,1928],
             # '06-01-20':[48,212,373,523,627],
             # '06-03-20':[44,77,84,154,313,463,731,822]
             # }
##### selected events with fast jumps #####
# whitelist = {'05-24-20':[193,289,422,478,1347,1510],
             # '05-26-20':[560],
             # '05-27-20':[411],
             # '05-28-20':[1445,1705],#1929
             # '05-29-20':[495,1514,1569], #325
             # '05-30-20':[34,787],
             # '05-31-20':[1033],
             # '06-01-20':[48,373,523]
             # }
##### 
whitelist = { '06-03-20':[1288,1343,1669],
              '06-04-20':[419, 569, 602, 612, 732, 780, 986, 987, 1025, 1074, 
                          1079,1098,  1374],
              '06-05-20':[11, 24, 55, 118, 204, 205, 226, 420, 690, 775, 847,
                          1048,  1192,], #1192
              '06-06-20':[92,287,289,464,511,523,950],
              '06-07-20':[428,461,753,1158,1727,1810,1921],
              '06-08-20':[22,  256,282,652,775,1052,1057,1302,1400,1478,1776,
                          1877,1980,2004],
              '06-09-20':[27,29,30,465,514,750,810,903,977],
              '07-11-20':[2166.7, 2858.0],
              '07-12-20':[51.5], # 515.1, 1590.1
              '07-13-20':[2156.6],
              '07-13-20':[154.1, 1344.5, 1557.2, 1587.5],
              '07-14-20':[2953.0], # 1990.3, 2346.7, 2933.5, 3903.2
              # '07-15-20':[95.3, 100.8, 664.0, 719.7], # 118.0,152.5,168.6,224.0,319.0
              # '07-16-20':[227.0, 349.0, 445.8, 518.9, 831.5, 871.8, 960.0, 
                          # 1008.9, 1200.5, 1210.1, 1292.5, 1833.1, 1857.5, 
                          # 1860.6, 1861.4, 1869.1], # 
              '07-15-20':[95.3, 100.8, 719.7], # 118.0,152.5,168.6,224.0,319.0
              '07-16-20':[871.8, 1200.5], # 
              '07-19-20':[32.0],
              '07-19-20':[],
              '07-21-20':[587.0, 1312.8, 1384.2, 1452.3, 2530.2],
              '07-22-20':[175.7, 7661.1],
             }
# weird blips: 7/12: 1275.6, 1815.4, 1892.9, 1899.2, 1901.9, 1904.0, 1916.8, 
                   # 2119.2, 2141.6
blacklist = {'05-29-20': { 13:[8], 544:[4], 576:[4], 1867:[1] },
             '05-30-20': { 34:[8,9], 787:[4] },
             '05-31-20': { 1033:[0], 1539:[2], 1928:[1,8] },
             '06-01-20': { 48:[5,6], 212:[2], 373:[6] },
             '06-03-20': { 44:[3,8], 77:[6], 84:[7], 463:[3,5,7], 822:[4], 1288:[0], 1343:[6] },
             '06-04-20': { 419:[1], 612:[7], 732:[1], 780:[3,9], 1025:[6], 1074:[4], 1374:[2,6,7] },
             '06-05-20': { 24:[4], 118:[4,6], 226:[3], 775:[1,4], 847:[1], 1192:[5,7] },
             '06-06-20': { 92:[1], 287:[2,5,7], 464:[7], 511:[4,6], 523:[9], 950:[7] },
             '06-07-20': { 428:[9], 461:[3], 753:[1], 1727:[1] },
             '06-08-20': { 22:[3], 775:[1,2,9], 1302:[2,8,9], 1400:[6,9], 1478:[1,2,9], 1776:[4,6], 1877:[8], 1980:[3], 2004:[5,7] },
             '06-09-20': { 27:[4], 810:[2,7], 903:[6,8], 977:[1,8] }
            }
data_files = {}

print len(files)
for k,f in enumerate(files): #
    path, num = noiselib.path_to_num(f)
    date = f.split('\\')[-4]
    if k == 0:
        base_trace = noiselib.loadmat(path.format(num+2))['Single_Shot_Occupation_RO1_SB1']
    if (whitelist is not None 
        and (date not in whitelist or num not in np.floor(whitelist[date]))):
        continue
    DF = TwoMeasDataFile(path.format(num), 
                         charge_var='Single_Shot_Occupations_RO1_SB1',
                                  #'Single_Shot_Occupations_RO2_SB1'], 
                         meas_var=['Single_Shot_Occupations_RO2_SB2',
                                   'Single_Shot_Occupations_RO1_SB2'])
    # DF.set_trigger_params(1000000., 300, 0.5)
    DF.set_trigger_params(1000000., 300, .25)
    trials = None
    if whitelist is not None and np.array(whitelist[date]).dtype == float:
        i = np.argwhere( np.floor(whitelist[date]) == num )
        v = np.array(whitelist[date])[i]
        trials = np.round((10*(v%1))).reshape(-1).astype(int)
    trigs = DF.get_triggers(trials)
    trigs = [(t,r) for t,r in trigs if date not in blacklist
                                    or num not in blacklist[date]
                                    or t not in blacklist[date][num]]
    # data_files[f] = DF
    print k, date, num, trigs
    for trial, rep in trigs:
        DF.add_base_ramsey_trace(base_trace)
        # DF.plot(trial)
        # plt.show()
        M2 = DF.get_P1_around_trig((trial,rep), meas_index=0)
        M1 = DF.get_P1_around_trig((trial,rep), meas_index=1)
        p10 = (M1 == 1) & (M2 == 0)
        p01 = (M1 == 0) & (M2 == 1)
        p11 = (M1 == 1) & (M2 == 1)
        p00 = (M1 == 0) & (M2 == 0)
        p1 = (M1 == 1)
        p0 = (M1 == 0)
        
        nreps = M2.size
        n_trace[N-rep:N+nreps-rep] += 1
        P10[N-rep:N+nreps-rep] += p10
        P01[N-rep:N+nreps-rep] += p01
        P00[N-rep:N+nreps-rep] += p00
        P11[N-rep:N+nreps-rep] += p11
        n0[N-rep:N+nreps-rep] += p0
        n1[N-rep:N+nreps-rep] += p1

# n0,n1=1.,1.
t = np.arange(-N,N+1)/2.
fig, ax = plt.subplots(1,1)
# ax.plot( t, P11+P10, label='P11 + P10' )
# ax.plot( t, P00+P10, label='P00 + P10' )
# ax.plot( t, 1.*P10/n1, label='Averaged P10' )
ax.plot( t, 1.*P01/n0, label='Averaged P01' )
# ax.plot( t, 1.*P00/n0, label='Averaged P00' )
# ax.plot( t, 1.*P11/n1, label='Averaged P11' )
# ax.plot( t, noiselib.movingmean(1.*P10/n1, 30), label='Averaged P10 smoothed 30' )
ax.plot( t, noiselib.movingmean(1.*P01/n0, 30), label='Averaged P01 smoothed 30' )
ax.plot( t, 1.*n_trace, 'k', label='number of files' )

ax.set_xlabel('Time From Trigger [ms]')
ax.legend()

plt.draw()
plt.pause(0.05)
