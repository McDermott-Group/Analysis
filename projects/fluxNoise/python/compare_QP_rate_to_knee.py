import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean
import QPTunneling
reload(QPTunneling)
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
files_charge = []
files_P1 = []
files_XIY_upper = []
files_XIY_lower = []
charge_path = ('fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General'
               '\Parameter\{}_correlation.hdf5')
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\{}\\{}\MATLABData\\'

# for d in [ds.q1q2_0701_charge_QP, ds.q1q2_0702_charge_QP]: # pre-calibration only
# for d in [ds.q1q2_0703_charge_QP, ds.q1q2_0704_charge_QP, ds.q1q2_0706_charge_QP]: # pre,post cal, 500us
# for d in [ds.q1q2_0707_charge_QP, ds.q1q2_0708_charge_QP]: # 30us
for d in [ds.q1q2_0710_charge_QP_2]: 
          
    CO.add_dataset(charge_path.format(d['path_charge']))
    if 'exclude_data' in d:
        for start,end in d['exclude_data']:
            CO.limit_dataset('Q1',start=start,end=end)
        
    files_charge += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                [path.format(date,'Charge_resetting') for date in d['date']], 1.)
    files_P1 += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                [path.format(date,'P1_I_X') for date in d['date']], 1.)
    files_XIY_upper += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                [path.format(date,'Fidelity_UpperBand') for date in d['date']], 1.)
    files_XIY_lower += CO.get_files_triggered_on_bad_fit(d['path_charge'], 'Q{}'.format(Q_A), 
                [path.format(date,'Fidelity_LowerBand') for date in d['date']], 1.)
  
  
P1_I = []
P1_X = []
XIY_upper = []
XIY_lower = []
flips0 = np.array([])
flips = np.array([])

print len(files)
for k,f in enumerate(files_charge[:50]):
    path, num = noiselib.path_to_num(f)
    # num = num + Q_A
    files_charge[k] = path.format(num+1)
    DF = TwoMeasDataFile(files_charge[k], charge_var='Single_Shot_Occupations_RO2_SB1', 
                                          meas_var=['Single_Shot_Occupations_RO2_SB2',
                                                    'Single_Shot_Occupations_RO1_SB2'],
                                          meas_fn=np.logical_xor)
    data = noiselib.loadmat(files_P1[k])
    P1_I += [data['Single_Shot_Occupation_RO1_SB2']]
    P1_X += [np.mean(np.logical_xor(data['Single_Shot_Occupations_RO1_SB2'],
                                   data['Single_Shot_Occupations_RO2_SB2']))]
    data = noiselib.loadmat(files_XIY_upper[k])
    XIY_upper += [np.mean(np.logical_xor(data['Single_Shot_Occupations_RO1_SB2'],
                                        data['Single_Shot_Occupations_RO2_SB2']))]
    data = noiselib.loadmat(files_XIY_lower[k])
    XIY_lower += [np.mean(np.logical_xor(data['Single_Shot_Occupations_RO1_SB2'],
                                        data['Single_Shot_Occupations_RO2_SB2']))]
    
    flips0 = np.append(flips0, np.sum( np.abs(np.diff(DF.o_meas,axis=1)), axis=1) )
    DF.apply_infidelity_correction_HMM(fidelity=[1.-XIY_lower[-1], XIY_upper[-1]])
    # DF.apply_infidelity_correction_HMM(fidelity=[0.95,0.95])
    flips = np.append(flips, np.sum( np.abs(np.diff(DF.o_meas,axis=1)), axis=1) )
    # DF.set_trigger_params(1000000., 300, 0.5)
    # DF.plot(5)

fig,(ax1,ax2,ax3) = plt.subplots(3, 1, sharex=True)
ax1.plot(P1_I, label='P1_I')
ax1.plot(P1_X, label='P1_X')
ax1.plot(XIY_upper, label='XIY_upper')
ax1.plot(XIY_lower, label='XIY_lower')
ax1.legend()
ax2.plot(np.arange(0,len(flips)/10.,0.1), flips0/8192/30e-6, label='flip rate')
ax2.plot(np.arange(0,len(flips)/10.,0.1), flips/8192/30e-6, label='corrected flip rate')
ax2.legend()
ax3.plot( np.abs(CO.get_jump_sizes()[0]['Q1']), label='jump size' )
ax3.plot( np.abs(CO.fit_R2[d['path_charge']]['Q1']), label='fit R^2' )
ax3.legend()
# ax.set_xlabel('File')
# ax.set_ylabel('P1')
plt.draw()
plt.pause(0.05)

QPT = QPTunneling(fs=1/30e-6)
QPT.add_datasets(files_charge[0:30], 'Single_Shot_Occupations_RO2_SB2')
QPT.plot_psd(figNum=111)