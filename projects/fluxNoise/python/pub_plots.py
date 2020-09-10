import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.special import erfc


# for publication
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 1.
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
fig_path = r'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs'

with open('dump_T1_sums.dat', 'rb') as f:
    P10,P01,P00,P11,n_trace,n0,n1 = pickle.load(f)
N = 10000
t = np.arange(-N,N+1)/2.
    
"""Plot T1 dropouts in occupation."""
# def gaus(x, b, A, sigma):
    # mu = 0
    # return b + A * exp(-(x-mu)**2/(2*sigma**2))
# def gausexp(x, b, A, sigma, A2, tau):
    # mu = 0
    # pulse = np.zeros(x.size)
    # pulse[x.size/2:] += A2*exp(-x[x.size/2:]/tau)
    # return gaus(x, b, A, sigma) + pulse
# def gausexp2(x, b, A, sigma, tau):
    # return b + 0.5 * A * exp((sigma**2-2*x*tau)/(2*tau**2)) * erfc((sigma**2-x*tau)/(np.sqrt(2)*sigma*tau))
# def gausexp3(x, b1, b2, A1, A2, sigma, tau1, tau2):
    # return np.concatenate([ gausexp2(x[:x.size/2],b1,A1,sigma,tau1),
                            # gausexp2(x[x.size/2:],b2,A2,sigma,tau2) ])
# fig = plt.figure(figsize=(3.375,2))
# # ax = fig.add_axes([0.1,0.1,1,1])
# ax = fig.add_subplot(111)
# ax.set_xlim(-100,100)
# ax.set_ylim(0.4,0.75)
# ax.set_xlabel('Time From Trigger (ms)')
# ax.set_ylabel('Average Occupation')
# for q in ['Q2','Q4']:
    # o = noiselib.movingmean(1.*P01[q]/n0[q], 1)
    # # popt, pcov = curve_fit(gaus, t[8000:12000], o[8000:12000], p0=[0.6,o[8000:12000].min()-0.6,5.])
    # # popt, pcov = curve_fit(gausexp, t[8000:12000], o[8000:12000], 
                            # # p0=[0.6,o[8000:12000].min()-0.6,5.,-0.3,5])
    # popt, pcov = curve_fit(gausexp2, t[9000:11000], o[9000:11000], 
                            # p0=[0.6,-0.1,2.,12.])
    # # popt, pcov = curve_fit(gausexp3, t[9000:11000], o[9000:11000], 
                            # # p0=[0.6,-0.1,2.,12.])
    # print popt
    # ax.plot( t, o, label='Averaged P01 smoothed 5 {}'.format(q) )
    # ax.plot( t, gausexp2(t,*popt), label='Averaged P01 smoothed 5 {}'.format(q) )
# fig.tight_layout()
# fig.savefig(fig_path+'\T1_occupation.pdf')

"""Plot T1 dropouts in T1 time."""
# fig = plt.figure(figsize=(3.375,2))
# # ax = fig.add_axes([0.1,0.1,1,1])
# ax = fig.add_subplot(111)
# ax.set_xlim(-100,100)
# ax.set_ylim(0,80)
# ax.set_xlabel('Time From Trigger (ms)')
# ax.set_ylabel('Inferred T1 Time ($\mu$s)')
# # a = {'Q2':0.6551, 'Q4':0.708}
# # c = {'Q2':0.1785, 'Q4':0.06842}
# a = {'Q2':0.5385, 'Q4':0.5166}
# c = {'Q2':0.2194, 'Q4':0.232}
# for q in ['Q2','Q4']:
    # o = noiselib.movingmean(1.*P01[q]/n0[q], 5)
    # T1 = -10./np.log((o-c[q])/a[q])
    # ax.plot( t, T1, label='Averaged P01 smoothed 5 {}'.format(q) )
# fig.tight_layout()
# fig.savefig(fig_path+'\T1_T1.pdf')

"""Plot single charge jump trace"""
# # 2166.7, 1069.3, 1830.0, 1256.5, *497.8*, 73.8, 2104.8
# # for k,f in enumerate(files[:]): # 290
    # # path, num = noiselib.path_to_num(f)
    # # if num == 497:
        # # print f
# p = (r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q4Corr\General\08-04-20\Charge_resetting\MATLABData\Charge_resetting_497.mat')
# data = noiselib.loadmat(p)
# charge_trace = data['Single_Shot_Occupations_RO2_SB1']
# fig = plt.figure(figsize=(3.375,2))
# # ax = fig.add_axes([0.1,0.1,1,1])
# ax = fig.add_subplot(111)
# # ax.set_xlim(-200,200)
# ax.set_xlabel('Time From Trigger')
# ax.set_ylabel('Average Occupation\nOver 30 Measurements')
# trig = 1570
# t = 1 * np.arange(8000)-trig
# ax.plot(t, noiselib.movingmean(charge_trace[8], 30) )
# ax.arrow(0,1,0,-0.2, width=20, head_width=100, head_length=0.05, fc='black',
        # length_includes_head=True)
# ax.plot([0],[noiselib.movingmean(charge_trace[8], 30)[trig]], 'xr', markersize=5)
# fig.tight_layout()
# fig.savefig(fig_path+'\charge_jump.pdf')

"""Plot Charge Offset Large"""
charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
CO = ChargeOffset()
for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
    CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
# CO.plot_charge_offset()
time, offset = CO.get_charge_offset()
fig, ax = plt.subplots(1, 1, figsize=(5.25,2.5), constrained_layout=True)
t0 = 200
time = (time - time[t0])/60./60.
for q in ['Q1','Q2','Q3','Q4']:
    offset[q] = offset[q] - offset[q][t0]
    ax.plot( time, offset[q], label=q)#, linewidth=1. )
ax.set_xlim(time[t0], time[t0]+10)
trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
y_max = np.vstack(offset.values())[:,trange].max()
y_min = np.vstack(offset.values())[:,trange].min()
ax.set_ylim(y_min, y_max)
rect = patches.Rectangle((1,1.5),0.75,1.5,linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)
# ax.legend()
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Charge Offset (e)')
# fig.tight_layout()
fig.savefig(fig_path+'\offset_charge.pdf')

"""Plot Charge Offset Zoom"""
charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
CO = ChargeOffset()
for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
    CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
# CO.plot_charge_offset()
time, offset = CO.get_charge_offset()
fig, ax = plt.subplots(1, 1, figsize=(1.95,2.5), constrained_layout=True)
t0 = 200
time = (time - time[t0])/60.
for q in ['Q1','Q2','Q3','Q4']:
    offset[q] = offset[q] - offset[q][t0]
    ax.plot( time, offset[q], label=q)#, linewidth=1. )
ax.set_xlim(60, 105)
trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
y_max = np.vstack(offset.values())[:,trange].max()
y_min = np.vstack(offset.values())[:,trange].min()
ax.set_ylim(1.5, 3)
# ax.legend()
ax.set_xlabel('Time (minutes)')
# ax.set_ylabel('Charge Offset (e)')
# fig.tight_layout()
fig.savefig(fig_path+'\offset_charge_zoom.pdf')



plt.draw()
plt.pause(0.05)