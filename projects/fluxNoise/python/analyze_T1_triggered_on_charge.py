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
import re
import pickle

# Q_A, Q_B = 1,2
# CO = ChargeOffset()
# files = []
# charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\{}Corr\General'
               # '\Parameter\{}_correlation.hdf5')
# path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\{}Corr\\General\\{}\\{}\MATLABData\\'

# # for d in [ds.q1q2_0524_charge_T1, ds.q1q2_0524_charge_T1_2, ds.q1q2_0525_charge_T1, 
          # # ds.q1q2_0526_charge_T1, ds.q1q2_0526_charge_T1_2, ds.q1q2_0527_charge_T1, 
          # # ds.q1q2_0528_charge_T1, ds.q1q2_0528_charge_T1_2]: # 200us cycle, 3us wait
# # for d in [ds.q1q2_0528_charge_T1_3, ds.q1q2_0529_charge_T1, ds.q1q2_0530_charge_T1, 
          # # ds.q1q2_0603_charge_T1]: # 500us cycle, 3us wait
# # for d in [ds.q1q2_0603_charge_T1_2, ds.q1q2_0604_charge_T1, ds.q1q2_0605_charge_T1,
          # # ds.q1q2_0607_charge_T1, ds.q1q2_0608_charge_T1]: # 500us cycle, 9us wait
# # for d in [ds.q1q2_0711_charge_T1, ds.q1q2_0714_charge_T1, ds.q1q2_0715_charge_T1,
            # # ds.q1q2_0719_charge_T1, ds.q1q2_0720_charge_T1, ds.q1q2_0722_charge_T1,
            # # ds.q1q2_0722_charge_T1_2, ds.q1q2_0724_charge_T1, ds.q1q2_0724_charge_T1_2,
            # # ds.q1q2_0725_charge_T1]: # 30us duty cycle, 10us wait, Q1-Q2
# # for d in [ds.q1q2q4_0726_charge_T1, ds.q1q2q4_0727_charge_T1, ds.q1q2q4_0728_charge_T1,
          # # ds.q1q2q4_0729_charge_T1, ds.q1q2q4_0730_charge_T1, ds.q1q2q4_0731_charge_T1, 
          # # ds.q1q2q4_0801_charge_T1, ds.q1q2q4_0802_charge_T1, ds.q1q2q4_0803_charge_T1]: # 30us duty cycle, 10us wait, Q1-Q2-Q4
# # for d in [ds.q1q4_0803_charge_T1, ds.q1q4_0805_charge_T1]: # 30us duty cycle, 10us wait, Q1-Q4 

# for d in [ds.q1q2_0711_charge_T1, ds.q1q2_0714_charge_T1, ds.q1q2_0715_charge_T1,
            # ds.q1q2_0719_charge_T1, ds.q1q2_0720_charge_T1, ds.q1q2_0722_charge_T1,
            # ds.q1q2_0722_charge_T1_2, ds.q1q2_0724_charge_T1, ds.q1q2_0724_charge_T1_2,
            # ds.q1q2_0725_charge_T1, ds.q1q2q4_0726_charge_T1, ds.q1q2q4_0727_charge_T1, 
            # ds.q1q2q4_0728_charge_T1, ds.q1q2q4_0729_charge_T1, ds.q1q2q4_0730_charge_T1, 
            # ds.q1q2q4_0731_charge_T1, ds.q1q2q4_0801_charge_T1, ds.q1q2q4_0802_charge_T1, 
            # ds.q1q4_0803_charge_T1, ds.q1q4_0805_charge_T1,
            # ds.q1q2q4_0810_charge_T1, ds.q1q2q4_0811_charge_T1, ds.q1q2q4_0812_charge_T1]: # ALL
# # for d in [ds.q1q2q4_0812_charge_T1]: #new 
          
    # CO.add_dataset(charge_path.format(d['Q'], d['path_charge']))
    # if 'exclude_data' in d:
        # for start,end in d['exclude_data']:
            # CO.limit_dataset('Q1',start=start,end=end)
        
    # files += CO.get_files_triggered_on_bad_fit( d['path_charge'], 'Q{}'.format(Q_A), 
                # [path.format(d['Q'], date,'Charge_resetting') for date in d['date']], 
                # d['thresh_charge'])
# with open('dump_T1_files.dat', 'wb') as f:
    # pickle.dump(files, f)

    
N = 10000
P10 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
P01 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
P00 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
P11 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
n_trace = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
n0 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
n1 = {'Q2': np.zeros(2*N+1), 'Q4': np.zeros(2*N+1)}
M1_before_trig = {'Q2': np.array([]), 'Q4': np.array([])}


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
              '07-12-20':[51.5, 775.0, 1631.4, 1978.6], # 515.1, 1590.1
              '07-13-20':[164.5, 590.9, 614.1, 905.1, 934.1, 942.8, 1771.3, 2156.6, 2878.1],
              '07-14-20':[47.6, 2397.4, 2402.1, 2421.6, 2953.0, 3318.5, 
                          3680.6, 3903.2], # 1990.3, 2346.7, 2933.5, 3903.2
              '07-15-20':[95.3, 100.8, 169.6, 232.8, 263.5, 664.0, 719.7], # 118.0,152.5,168.6,224.0,319.0
              '07-16-20':[28.6, 59.2, 100.1, 158.7, 221.9, 349.0, 871.8, 1200.5], 
              '07-19-20':[32.0, 820.5, 846.4, 1183.1],
              '07-20-20':[538.4, 556.4, 910.8, 1352.3, 1463.5, 1528.1, 1681.8, 
                          2330.6, 2598.4, 2801.6],
              '07-21-20':[236.7, 493.7, 587.0, 815.0, 1312.8, 1384.2, 1452.3, 
                          1625.4, 1756.2, 1773.5, 2095.9, 2263.9, 2530.2, 
                          2732.6, 2788.9, 2910.4],
              '07-22-20':[175.7, 253.5, 420.8, 720.4, 739.4, 2316.4, 7661.1], # 1908.0
              '07-23-20':[271.6, 1069.3, 1248.9, 1334.3, 1467.1, 1467.6, 1830.0, 
                          1831.8, 1914.2, 1984.4, 2157.9, 2559.6, 2606.9, 2681.9],
              '07-24-20':[62.4, 107.5, 128.4, 367.5, 605.4, 614.4, 627.6, 727.0, 
                          968.0, 973.9, 1154.1, 1271.5, 1315.5, 1354.5, 1753.7, 
                          1762.6, 1914.7, 2063.3, 2352.5, 2385.5, 2551.4, 2602.4],
              '07-25-20':[18.1, 96.6, 321.3, 507.8, 545.7, 743.7, 1041.5, 1132.4,
                          1253.2, 1256.5, 1319.5, 1495.1, 1516.1, 2689.2],
              '07-26-20':[180.9, 505.1, 643.2, 662.4, 723.5],
              '07-27-20':[978.9, 2098.1],
              '07-28-20':[296.4, 314.8, 363.6, 786.7, 865.8, 907.5, 1091.3, 
                          1711.9, 1775.4, 1898.1],
              '07-29-20':[53.4, 141.0, 174.4, 702.0, 752.2, 971.3, 997.5, 1067.0,
                          1070.5, 1073.9, 1356.3, 1366.0, 1505.6, 1611.1, 1798.1,
                          1989.1, 2155.0],
              '07-30-20':[587.1, 938.2, 994.5, 1479.4, 1544.6],
              '07-31-20':[123.4, 189.6, 724.0, 855.4, 859.1, 1013.5, 1062.8, 
                          1391.6, 1660.6, 1662.8, 1906.6, 2274.4], # 1753.4, 
              '08-01-20':[44.5, 48.8, 67.1, 231.8, 289.8, 486.6, 624.7, 766.4,
                          774.6, 782.1, 867.1, 959.5, 980.3, 1126.2], # 1079.6
              '08-02-20':[70.1, 300.3, 442.1, 620.1, 817.3, 928.1, 1059.0, 
                          1093.1, 1117.0, 1269.5, 1429.8, 1938.4], 
              '08-03-20':[84.1, 110.7, 208.9, 320.6, 350.4, 363.6, 375.1, 460.1, 
                          552.7, 574.4, 692.4, 1141.6, 1184.7], 
              '08-04-20':[86.4, 123.7, 183.1, 254.7, 396.6, 497.8, 532.7, 
                          692.4, 772.6, 1099.5, 1333.0, 1479.6, 1652.1, 1792.4,
                          2089.1, 2121.4, 2168.5, 2246.8, 2252.1, 2277.4, 2351.6,
                          2432.2, 2451.8, 2497.7, 2735.9], 
              '08-05-20':[15.3, 73.8, 249.1, 255.5, 629.6, 669.66, 780.9, 920.7,
                          1010.2, 1112.5, 1117.7, 1153.1, 1175.7, 1511.1, 1526.6,
                          1901.4, 2104.8, 2106.8, 2167.1, 2249.5, 2883.5],
              '08-06-20':[33.6, 201.4, 452.6, 462.7, 541.9, 605.7, 628.4, 1009.4, 
                          1151.4, 1182.2, 1308.6, 1375.8, 1505.5, 1686.4], 
              # '08-06-20':[356.6, 366.6, 721.5, 733.2, 746.4, 760.9, 763.8, 948.1], 
              # '08-07-20':[185.7, 266.1, 361.5, 495.3, 849.3, 942.0, 1012.0, 1055.7,
                          # 1093.0, 1292.3, 1343.8, 1354.4, 1403.0, 1858.5, 1893.7,
                          # 2317.5, 2350.0, 2512.1, 2535.4, 2535.9, 2570.4], 
              '08-10-20':[601.7, 662.1, 704.6, 847.4, 1571.6, 1610.0],
              '08-11-20':[95.8, 129.6, 542.5],
              '08-12-20':[35.6, 440.6, 822.1],
              '08-13-20':[70.3, 631.4, 651.4, 880.9, 1016.2],
             }
# whitelist = None
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
reg = re.compile('Q1(Q\d)(Q\d)*')


with open('dump_T1_files.dat', 'rb') as f:
    files = pickle.load(f)
print len(files)
for k,f in enumerate(files[:]): # 290
    Qs = list(reg.findall(f)[0])
    if Qs[-1] == '':
        Qs.pop()
    path, num = noiselib.path_to_num(f)
    date = f.split('\\')[-4]
    if k == 0:
        base_trace = noiselib.loadmat(path.format(num+2))['Single_Shot_Occupation_RO2_SB1']
    if (whitelist is not None 
        and (date not in whitelist or num not in np.floor(whitelist[date]))):
        continue
    meas_vars = ['Single_Shot_Occupations_RO2_SB2',
                 'Single_Shot_Occupations_RO1_SB2']
    if len(Qs) > 1:
        meas_vars += ['Single_Shot_Occupations_RO2_SB3',
                      'Single_Shot_Occupations_RO1_SB3']
    DF = TwoMeasDataFile(path.format(num), 
                         charge_var='Single_Shot_Occupations_RO2_SB1',
                                  #'Single_Shot_Occupations_RO2_SB1'], 
                         # meas_var=['Single_Shot_Occupations_RO2_SB2',
                                   # 'Single_Shot_Occupations_RO1_SB2'])
                         # meas_var=['Single_Shot_Occupations_RO2_SB3',
                                   # 'Single_Shot_Occupations_RO1_SB3'])
                         meas_var=meas_vars)
    IsQsFile = noiselib.loadmat(path.format(num))
    # DF.set_trigger_params(1000000., 300, 0.5)
    DF.set_trigger_params(10000000., 500, 1)
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
        # plt.figure()
        # plt.scatter(IsQsFile['Is_RO2_SB2'][trial], IsQsFile['Qs_RO2_SB2'][trial], s=1)
        # plt.scatter(IsQsFile['Is_RO2_SB2'][trial][rep-10:rep+10], IsQsFile['Qs_RO2_SB2'][trial][rep-10:rep+10], c=cm.brg(np.arange(20)/20.))
        # plt.draw(); plt.pause(0.05);
        # DF.plot(trial)
        # plt.show()
        for i,q in enumerate(Qs):
            M2 = DF.get_P1_around_trig((trial,rep), meas_index=0+2*i)
            M1 = DF.get_P1_around_trig((trial,rep), meas_index=1+2*i)
            p10 = (M1 == 1) & (M2 == 0)
            p01 = (M1 == 0) & (M2 == 1)
            p11 = (M1 == 1) & (M2 == 1)
            p00 = (M1 == 0) & (M2 == 0)
            p1 = (M1 == 1)
            p0 = (M1 == 0)
            
            nreps = M2.size
            n_trace[q][N-rep:N+nreps-rep] += 1
            P10[q][N-rep:N+nreps-rep] += p10
            P01[q][N-rep:N+nreps-rep] += p01
            P00[q][N-rep:N+nreps-rep] += p00
            P11[q][N-rep:N+nreps-rep] += p11
            n0[q][N-rep:N+nreps-rep] += p0
            n1[q][N-rep:N+nreps-rep] += p1
            M1_before_trig[q] = np.concatenate( [M1_before_trig[q], M1[:rep]] )

with open('dump_T1_sums.dat', 'wb') as f:
    pickle.dump( (P10,P01,P00,P11,n_trace,n0,n1,M1_before_trig), f)



with open('dump_T1_sums.dat', 'rb') as f:
    P10,P01,P00,P11,n_trace,n0,n1,M1_before_trig = pickle.load(f)
# n0,n1=1.,1.
t = np.arange(-N,N+1)/2.
fig, ax = plt.subplots(1,1)
for q in ['Q2','Q4']:
    # ax.plot( t, P11+P10, label='P11 + P10' )
    # ax.plot( t, P00+P10, label='P00 + P10' )
    # ax.plot( t, 1.*P10/n1, label='Averaged P10' )
    ax.plot( t, 1.*P01[q]/n0[q], label='Averaged P01 {}'.format(q), linewidth=0.5 )
    # ax.plot( t, 1.*P00/n0, label='Averaged P00' )
    # ax.plot( t, 1.*P11/n1, label='Averaged P11' )
    # ax.plot( t, noiselib.movingmean(1.*P10/n1, 30), label='Averaged P10 smoothed 30' )
    ax.plot( t, noiselib.movingmean(1.*P01[q]/n0[q], 5), label='Averaged P01 smoothed 30 {}'.format(q) )
    ax.plot( t, 1.*n_trace[q], 'k', label='number of files {}'.format(q) )
ax.set_xlabel('Time From Trigger [ms]')
ax.legend()

plt.draw()
plt.pause(0.05)

print np.mean(M1_before_trig['Q2'])
print np.mean(M1_before_trig['Q4'])
