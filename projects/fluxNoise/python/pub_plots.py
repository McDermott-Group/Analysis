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


# for publication
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 10

with open('dump_T1_sums.dat', 'rb') as f:
    P10,P01,P00,P11,n_trace,n0,n1 = pickle.load(f)
N = 10000
t = np.arange(-N,N+1)/2.
    
"""Plot T1 dropouts in occupation."""
fig = plt.figure(figsize=(3.375,2))
# ax = fig.add_axes([0.1,0.1,1,1])
ax = fig.add_subplot(111)
ax.set_xlim(-100,100)
ax.set_ylim(0.4,0.75)
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.set_xlabel('Time From Trigger (ms)')
ax.set_ylabel('Average Occupation')
for q in ['Q2','Q4']:
    o = noiselib.movingmean(1.*P01[q]/n0[q], 5)
    ax.plot( t, o, label='Averaged P01 smoothed 5 {}'.format(q) )
fig.tight_layout()

"""Plot T1 dropouts in T1 time."""
fig = plt.figure(figsize=(3.375,2))
# ax = fig.add_axes([0.1,0.1,1,1])
ax = fig.add_subplot(111)
ax.set_xlim(-100,100)
ax.set_ylim(0,80)
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.set_xlabel('Time From Trigger (ms)')
ax.set_ylabel('Inferred T1 Time ($\mu$s)')
# a = {'Q2':0.6551, 'Q4':0.708}
# c = {'Q2':0.1785, 'Q4':0.06842}
a = {'Q2':0.5385, 'Q4':0.5166}
c = {'Q2':0.2194, 'Q4':0.232}
for q in ['Q2','Q4']:
    o = noiselib.movingmean(1.*P01[q]/n0[q], 5)
    T1 = -10./np.log((o-c[q])/a[q])
    ax.plot( t, T1, label='Averaged P01 smoothed 5 {}'.format(q) )
fig.tight_layout()

"""Plot single charge jump trace"""
# 2166.7, 1069.3, 1830.0, 1256.5, *497.8*, 73.8, 2104.8
# for k,f in enumerate(files[:]): # 290
    # path, num = noiselib.path_to_num(f)
    # if num == 497:
        # print f
p = (r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q4Corr\General\08-04-20\Charge_resetting\MATLABData\Charge_resetting_497.mat')
data = noiselib.loadmat(p)
charge_trace = data['Single_Shot_Occupations_RO2_SB1']
fig = plt.figure(figsize=(3.375,2))
# ax = fig.add_axes([0.1,0.1,1,1])
ax = fig.add_subplot(111)
# ax.set_xlim(-200,200)
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.set_xlabel('Time From Trigger')
ax.set_ylabel('Average Occupation\nOver 30 Measurements')
trig = 1570
t = 1 * np.arange(8000)-trig
ax.plot(t, noiselib.movingmean(charge_trace[8], 30) )
ax.arrow(0,1,0,-0.2, width=20, head_width=100, head_length=0.05, fc='black',
        length_includes_head=True)
ax.plot([0],[noiselib.movingmean(charge_trace[8], 30)[trig]], 'xr', markersize=5)
fig.tight_layout()

plt.draw()
plt.pause(0.05)
