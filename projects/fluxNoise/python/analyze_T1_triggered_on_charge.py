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
# CO.add_dataset('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
                # '\Parameter\cxs2305ngo_correlation.hdf5')
CO.add_dataset(charge_path.format('cxu0723thw')) # 4-5 events
CO.add_dataset(charge_path.format('cxu1901pfy'))
CO.add_dataset(charge_path.format('cxv0747mel'))
CO.add_dataset(charge_path.format('cxw1630xbt'))
CO.add_dataset(charge_path.format('cxx0332bgf')).limit_dataset('Q2',start=915)
CO.add_dataset(charge_path.format('cxx1757kfx'))
CO.add_dataset(charge_path.format('cxy0502hxc')).limit_dataset('Q2',start=905)
CO.add_dataset(charge_path.format('cxy2033hub'))
CO.add_dataset(charge_path.format('cxz0428jaq'))
CO.add_dataset(charge_path.format('cxz2301mxf'))
CO.add_dataset(charge_path.format('cyb0359hkw'))
CO.add_dataset(charge_path.format('cye0625zge'))
CO.add_dataset(charge_path.format('cye2324ghg'))
CO.add_dataset(charge_path.format('cyf2353zve'))
CO.add_dataset(charge_path.format('cyg2351lgx'))
CO.add_dataset(charge_path.format('cyi0501krn'))
CO.add_dataset(charge_path.format('cyj0526hzn'))

files = []
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\Charge_resetting\MATLABData\\'

######## 200us ########
# # files += CO.get_files_triggered_on_bad_fit('cxs2305ngo', 'Q{}'.format(Q_A), [path.format('05-22-20'), path.format('05-23-20')], thresh=0.8 )
# files += CO.get_files_triggered_on_bad_fit('cxu0723thw', 'Q{}'.format(Q_A), path.format('05-24-20'), thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxu1901pfy', 'Q{}'.format(Q_A), path.format('05-24-20'), thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxv0747mel', 'Q{}'.format(Q_A), path.format('05-25-20'), thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxw1630xbt', 'Q{}'.format(Q_A), path.format('05-26-20'), thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxx0332bgf', 'Q{}'.format(Q_A), [path.format('05-26-20'), path.format('05-27-20')], thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxx1757kfx', 'Q{}'.format(Q_A), path.format('05-27-20'), thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxy0502hxc', 'Q{}'.format(Q_A), path.format('05-28-20'), thresh=0.4 ) # pretty shitty dataset, blobs looked wrong
# files += CO.get_files_triggered_on_bad_fit('cxy2033hub', 'Q{}'.format(Q_A), path.format('05-28-20'), thresh=0.91 )

######## 500us ########
# files += CO.get_files_triggered_on_bad_fit('cxz0428jaq', 'Q{}'.format(Q_A), [path.format('05-28-20'),path.format('05-29-20')], thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cxz2301mxf', 'Q{}'.format(Q_A), [path.format('05-29-20'),path.format('05-30-20')], thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cyb0359hkw', 'Q{}'.format(Q_A), [path.format('05-30-20'),path.format('05-31-20'),path.format('06-01-20')], thresh=0.91 )
# files += CO.get_files_triggered_on_bad_fit('cye0625zge', 'Q{}'.format(Q_A), [path.format('06-03-20')], thresh=0.91 )

######## 10us idle ########
files += CO.get_files_triggered_on_bad_fit('cye2324ghg', 'Q{}'.format(Q_A), [path.format('06-03-20'),path.format('06-04-20')], thresh=0.91 )
files += CO.get_files_triggered_on_bad_fit('cyf2353zve', 'Q{}'.format(Q_A), [path.format('06-04-20'),path.format('06-05-20')], thresh=0.91 )
files += CO.get_files_triggered_on_bad_fit('cyg2351lgx', 'Q{}'.format(Q_A), [path.format('06-05-20'),path.format('06-06-20')], thresh=0.91 )
files += CO.get_files_triggered_on_bad_fit('cyi0501krn', 'Q{}'.format(Q_A), [path.format('06-07-20'),path.format('06-08-20')], thresh=0.91 )
files += CO.get_files_triggered_on_bad_fit('cyj0526hzn', 'Q{}'.format(Q_A), [path.format('06-08-20'),path.format('06-09-20')], thresh=0.91 )
    
P10 = np.zeros(2*3000+1)
P01 = np.zeros(2*3000+1)
P00 = np.zeros(2*3000+1)
P11 = np.zeros(2*3000+1)
n_trace = np.zeros(2*3000+1)
n0 = np.zeros(2*3000+1)
n1 = np.zeros(2*3000+1)


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
              '06-04-20':[419,569,602,612,732,780,986,987,1025,1074,1079,1098,  1374],
              '06-05-20':[11,24,55,118,204,205,226,420,690,775,847,1048,  1192,], #1192
              '06-06-20':[92,287,289,464,511,523,950],
              '06-07-20':[428,461,753,1158,1727,1810,1921],
              '06-08-20':[22,  256,282,652,775,1052,1057,1302,1400,1478,1776,1877,1980,2004],
              '06-09-20':[27,29,30,465,514,750,810,903,977]
             }
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
for k,f in enumerate(files): #400
    path, num = noiselib.path_to_num(f)
    date = f.split('\\')[-4]
    if (whitelist is not None 
        and (date not in whitelist or num not in whitelist[date])):
        continue
    DF = TwoMeasDataFile(path.format(num+Q_A-1), charge_var='Single_Shot_Occupations_RO2_SB1', 
                                                 meas_var=['Single_Shot_Occupations_RO2_SB2',
                                                           'Single_Shot_Occupations_RO1_SB2'])
    DF.set_trigger_params(1000000., 300, 0.5)
    # DF.set_trigger_params(1000000., 200, 1.1)
    # DF.apply_infidelity_correction(7)
    trigs = DF.get_triggers()
    trigs = [(t,r) for t,r in trigs if date not in blacklist
                                    or num not in blacklist[date]
                                    or t not in blacklist[date][num]]
    data_files[f] = DF
    print k, date, num, trigs
    for trial, rep in trigs:
        # DF.plot(trial)
        M2 = DF.get_P1_around_trig((trial,rep), meas_index=0)
        M1 = DF.get_P1_around_trig((trial,rep), meas_index=1)
        p10 = (M1 == 1) & (M2 == 0)
        p01 = (M1 == 0) & (M2 == 1)
        p11 = (M1 == 1) & (M2 == 1)
        p00 = (M1 == 0) & (M2 == 0)
        p1 = (M1 == 1)
        p0 = (M1 == 0)
        
        nreps = M2.size
        n_trace[3000-rep:3000+nreps-rep] += 1
        P10[3000-rep:3000+nreps-rep] += p10
        P01[3000-rep:3000+nreps-rep] += p01
        P00[3000-rep:3000+nreps-rep] += p00
        P11[3000-rep:3000+nreps-rep] += p11
        n0[3000-rep:3000+nreps-rep] += p0
        n1[3000-rep:3000+nreps-rep] += p1

# n0,n1=1.,1.
t = np.arange(-3000,3000+1)/2.
fig, ax = plt.subplots(1,1)
ax.plot( t, P11+P10, label='P11 + P10' )
ax.plot( t, P00+P10, label='P00 + P10' )
# ax.plot( t, 1.*P10/n1, label='Averaged P10' )
ax.plot( t, 1.*P01/n0, label='Averaged P01' )
ax.plot( t, 1.*P00/n0, label='Averaged P00' )
# ax.plot( t, 1.*P11/n1, label='Averaged P11' )
# ax.plot( t, noiselib.movingmean(1.*P10/n1, 30), label='Averaged P10 smoothed 30' )
# ax.plot( t, noiselib.movingmean(1.*P01/n0, 30), label='Averaged P01 smoothed 30' )
ax.plot( t, 1.*n_trace, 'k', label='number of files' )

ax.set_xlabel('Time From Trigger [ms]')
ax.legend()

plt.draw()
plt.pause(0.05)
