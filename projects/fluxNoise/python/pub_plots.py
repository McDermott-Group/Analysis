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
halfwidth = 3.5
fullwidth = 7.2

with open('dump_T1_sums.dat', 'rb') as f:
    P10,P01,P00,P11,n_trace,n0,n1,M1_before_trig = pickle.load(f)
N = 10000
t = 0.04*np.arange(-N,N+1)
    
"""Plot T1 dropouts in occupation."""
def gaus(x, b, A, sigma, mu=0):
    return b + A * exp(-(x-mu)**2/(2*sigma**2))
def gausexp(x, b, A, sigma, A2, tau):
    mu = 0
    pulse = np.zeros(x.size)
    pulse[x.size/2:] += A2*exp(-x[x.size/2:]/tau)
    return gaus(x, b, A, sigma) + pulse
def gausexp2(x, b, A, sigma, tau, mu=0):
    return b + 0.5 * A * np.exp((sigma**2-2*x*tau+2*mu*tau)/(2*tau**2)) * erfc((sigma**2-x*tau+mu*tau)/(np.sqrt(2)*sigma*tau))
def gausexp3(x, b1, b2, A1, A2, sigma, tau1, tau2):
    return np.concatenate([ gausexp2(x[:x.size/2],b1,A1,sigma,tau1),
                            gausexp2(x[x.size/2:],b2,A2,sigma,tau2) ])

fig = plt.figure(figsize=(halfwidth,2))
# ax = fig.add_axes([0.1,0.1,1,1])
ax = fig.add_subplot(111)
ax.set_xlim(-5,5)
ax.set_ylim(0.4,0.75)
# ax.set_xlabel('Time From Trigger (ms)')
ax.set_ylabel('$P_{01}$')
for i,q in enumerate(['Q2','Q4']):
    o = noiselib.movingmean(1.*P01[q]/n0[q], 1)
    # popt, pcov = curve_fit(gaus, t[8000:12000], o[8000:12000], p0=[0.6,o[8000:12000].min()-0.6,5.])
    # popt, pcov = curve_fit(gausexp, t[8000:12000], o[8000:12000], 
                            # p0=[0.6,o[8000:12000].min()-0.6,5.,-0.3,5])
    popt, pcov = curve_fit(gausexp2, t[9000:11000], o[9000:11000], 
                            p0=[0.6,-0.9,.3,.3],
                            bounds=(-np.inf,
                                    [np.inf,np.inf,np.inf,np.inf]))
                            # p0=[0.6,-0.1,2.,12.])
    # o2 = noiselib.movingmean(1.*P01['Q2']/n0['Q2'], 1)
    # o4 = noiselib.movingmean(1.*P01['Q4']/n0['Q4'], 1)
    # popt, pcov = curve_fit(gausexp3, np.concatenate([t[9000:11000],t[9000:11000]]), 
                            # np.concatenate([o2[9000:11000], o4[9000:11000]]), 
                            # p0=[0.6,0.6,-1,-1,2.,1.,1.])
    print popt
    ax.plot( t, o, label='Averaged P01 smoothed 5 {}'.format(q), color='C{}'.format(i+1) )
    ax.plot( t, gausexp2(t,*popt), label='Averaged P01 smoothed 5 {}'.format(q) )
    # ax.plot( np.concatenate([t,t]), gausexp3(np.concatenate([t,t]),*popt), label='Averaged P01 smoothed 5 {}'.format(q) )
fig.tight_layout()
fig.savefig(fig_path+'\T1_occupation.pdf')

"""Plot T1 dropouts in T1 time."""
fig = plt.figure(figsize=(halfwidth,2))
# ax = fig.add_axes([0.1,0.1,1,1])
ax = fig.add_subplot(111)
axr = ax.twinx()
ax.set_xlim(-5,5)
# ax.set_ylim(0,80)
ax.set_xlabel('Time From Trigger (ms)')
# ax.set_ylabel('Inferred T1 Time ($\mu$s)')
ax.set_ylabel('$\Delta\Gamma_{01}$ ($\mu s^{-1}$)')
axr.set_ylabel('$x_{QP}$')
# a = {'Q2':0.6551, 'Q4':0.708}
# c = {'Q2':0.1785, 'Q4':0.06842}
# a = {'Q2':0.5385, 'Q4':0.5166}
# c = {'Q2':0.2194, 'Q4':0.232}
# a = {'Q2':0.7592, 'Q4':0.5423}
# c = {'Q2':0.1439, 'Q4':0.1198}
a = {'Q2':1-2*0.38, 'Q4':1-2*0.31} # amp
c = {'Q2':0.38, 'Q4':0.31} # offset
df01 = {'Q2':(4.436-4.4308)/2, 'Q4':(4.4053-4.3995)/2}
for i,q in enumerate(['Q2','Q4']):
    o = noiselib.movingmean(1.*P01[q]/n0[q], 1)
    T1 = -10./np.log((o-c[q])/a[q])
    gamma = 1./T1
    x_qp = np.pi*gamma/np.sqrt( 2*2*np.pi*df01[q] )
    ax.plot( t, 1./T1, label='$P_{}$ {}'.format('01',q), color='C{}'.format(i+1) )
    axr.plot(t, x_qp)
fig.tight_layout()
fig.savefig(fig_path+'\T1_T1.pdf')

"""Plot single charge jump trace"""
# # 2166.7, 1069.3, 1830.0, 1256.5, *497.8*, 73.8, 2104.8
# # for k,f in enumerate(files[:]): # 290
    # # path, num = noiselib.path_to_num(f)
    # # if num == 497:
        # # print f
# p = (r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q4Corr\General\08-04-20\Charge_resetting\MATLABData\Charge_resetting_497.mat')
# data = noiselib.loadmat(p)
# charge_trace = data['Single_Shot_Occupations_RO2_SB1']
# fig = plt.figure(figsize=(halfwidth,2))
# # ax = fig.add_axes([0.1,0.1,1,1])
# ax = fig.add_subplot(111)
# # ax.set_xlim(-200,200)
# # ax.set_xlabel('Time From Trigger')
# ax.set_ylabel('$P_{01}$')
# trig = 1570
# t = 0.040 * (np.arange(8000)-trig)
# ax.plot(t, noiselib.movingmean(charge_trace[8], 30) )
# ax.arrow(0,1,0,-0.2, width=0.1, head_width=0.300, head_length=0.05, fc='black',
        # length_includes_head=True)
# ax.plot([0],[noiselib.movingmean(charge_trace[8], 30)[trig]], 'xr', markersize=5)
# ax.set_xlim(-5,5)
# fig.tight_layout()
# fig.savefig(fig_path+'\charge_jump.pdf')

"""Plot Charge Offset Large"""
# charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
# CO = ChargeOffset()
# for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
    # CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
# # CO.plot_charge_offset()
# time, offset = CO.get_charge_offset()
# fig, ax = plt.subplots(1, 1, figsize=(5.25,2.5), constrained_layout=True)
# t0 = 200
# time = (time - time[t0])/60./60.
# for q in ['Q1','Q2','Q3','Q4']:
    # offset[q] = offset[q] - offset[q][t0]
    # ax.plot( time, offset[q], label=q)#, linewidth=1. )
# ax.set_xlim(time[t0], time[t0]+10)
# trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
# y_max = np.vstack(offset.values())[:,trange].max()
# y_min = np.vstack(offset.values())[:,trange].min()
# ax.set_ylim(y_min, y_max)
# rect = patches.Rectangle((1,1.5),0.75,1.5,linewidth=1,edgecolor='k',facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((175./60,-1),20./60,2,linewidth=1,edgecolor='k',facecolor='none')
# ax.add_patch(rect)
# # ax.legend()
# ax.set_xlabel('Time (hours)')
# ax.set_ylabel('Charge Offset (e)')
# # fig.tight_layout()
# fig.savefig(fig_path+'\offset_charge.pdf')

"""Plot Charge Offset Zoom"""
# charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
# CO = ChargeOffset()
# for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
    # CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
# # CO.plot_charge_offset()
# time, offset = CO.get_charge_offset()
# fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(1.95,2.5), constrained_layout=True)
# t0 = 200
# time = (time - time[t0])/60.
# for q in ['Q1','Q2','Q3','Q4']:
    # offset[q] = offset[q] - offset[q][t0]
    # ax1.plot( time, offset[q], label=q)#, linewidth=1. )
    # ax2.plot( time, offset[q], label=q)#, linewidth=1. )
# ax1.set_xlim(60, 105)
# ax2.set_xlim(175, 195)
# trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
# y_max = np.vstack(offset.values())[:,trange].max()
# y_min = np.vstack(offset.values())[:,trange].min()
# ax1.set_ylim(1.5, 3)
# ax2.set_ylim(-1, 1)
# # ax.legend()
# ax2.set_xlabel('Time (minutes)')
# # ax.set_ylabel('Charge Offset (e)')
# # fig.tight_layout()
# fig.savefig(fig_path+'\offset_charge_zoom.pdf')

"""Plot Bloch Sphere"""



plt.draw()
plt.pause(0.05)