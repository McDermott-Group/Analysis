import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from noiselib import matpaths
from QPTunneling import *
import ChargeOffset
importlib.reload(ChargeOffset)
from ChargeOffset import *
from datasets import *

def get_averaged_P1(files, label=''):
    P1_avg = np.array([])
    for f in files:
        try:
            data = noiselib.loadmat(f)
        except:
            print('corrupted:', f)
            data = { 'Single_Shot_Occupation{}'.format(label): np.nan }
        o = np.array(data['Single_Shot_Occupation{}'.format(label)])
        P1_avg = np.append(P1_avg, o)
    return P1_avg
    
def get_averaged_flips(files, fidelity, label='',axis=None):
    flip_avg = np.array([])
    for i,f in enumerate(files):
        data = noiselib.loadmat(f)
        o = np.array(data['Single_Shot_Occupations{}'.format(label)])
        # o = noiselib.apply_infidelity_correction(o, 9)
        try:
            o = noiselib.apply_infidelity_correction_HMM(o, fidelity=fidelity[i])
        except:
            print('failed to correct')
            print(f)
            print(fidelity[i])
        flip_avg = np.append(flip_avg, np.mean(np.abs(np.diff(o, axis=1)),axis=axis))
    return flip_avg

def plot_data(data, xlabel='', ylabel='', title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.plot(data)
    plt.draw()
    plt.pause(0.05)
    return ax

####################################
###fluxNoise
# date, Q = '03-17-20', 'Q1'
# P1_files = (0,1003+1)
# QP_files = (0,1003+1)
# filter_levels = (0.05, 0.15) # low, high
# date, Q = '03-17-20', 'Q2'
# P1_files = (0,908+1)
# QP_files = (0,908+1)
# filter_levels = (0.05, 0.05) # low, high
# date, Q = '03-18-20', 'Q3'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.02, 0.02) # low, high
# date, Q = '03-18-20', 'Q4'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.15, 0.2) # low, high
### fluxNoise2
# date, Q = '03-17-20', 'Q1'
# P1_files = (0,1003+1)
# QP_files = (0,1003+1)
# filter_levels = (0.05, 0.15) # low, high
# date, Q = '04-10-20', 'Q2'
# P1_files = (0,999+1)
# QP_files = (0,999+1)
# filter_levels = (0.05, 0.05) # low, high

###fluxNoise2

# date, Q = ['04-28-20','04-29-20'], 'Q3'
# filesP1, fName = [range(0,821+1),range(0,110+1)], ['P1_I','P1_X']
# filesQP = [range(0,821+1),range(0,110+1)]

# date, Q = '04-30-20', 'Q4'
# filesP1, fName = range(234,470+1), ['P1_I','P1_X']
# filesQP = range(1,237+1)
# fileCharge = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q4\General\Parameter\cww1911obk_charge_offset.hdf5'

# date, Q = '05-01-20', 'Q3'
# filesP1, fName = range(0,664+1), ['P1_I','P1_X']
# filesQP = range(0,663+1)
# fileCharge = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q3\General\Parameter\cwx0615mbj_charge_offset.hdf5'

date, Q = '05-04-20', 'Q1'
filesP1, fName = list(range(0,552+1)), ['P1_I','P1_X']
filesQP = list(range(0,552+1)) #552
fileCharge = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1\General\Parameter\cxa1458apn_charge_offset.hdf5'

#####################################

# ds = q3_0428_QP
ds = q4_0430_QP
# ds = q3_0501_QP
# ds = q1_0504_QP

fP1_I = matpaths(fileName='P1_I', fileNums='files_P1', **ds)
fP1_X = matpaths(fileName='P1_X', fileNums='files_P1', **ds)
fQP = matpaths(fileName='QP_Tunneling_PSD', fileNums='files_QP', **ds)

P1_avg_I = get_averaged_P1(fP1_I, '')
P1_avg_X = get_averaged_P1(fP1_X, '')
flip_avg = get_averaged_flips(fQP, 
                fidelity=np.transpose([1-P1_avg_I, P1_avg_X]), label='',axis=1)
flip_avg[flip_avg < 10./8192] = np.nan
CO = ChargeOffset()
CO.add_dataset(fileCharge)
charge_jumps, sigma = CO.get_jump_sizes()
charge_jumps = np.abs(list(charge_jumps.values())[0])

fig, axs = plt.subplots(4, sharex=True) #, gridspec_kw={'hspace':5})
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
labels = ['P1_I','P1_X','Avg Flip Rate','Charge Jumps']
for i,data in enumerate([ P1_avg_I, P1_avg_X, flip_avg, charge_jumps ]):
    axs[i].plot(data, color=colors[i], label=labels[i])
    axs[i].legend()
axs[-1].set_xlabel('File')
plt.draw()
plt.pause(0.05)