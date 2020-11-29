import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
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
from scipy import asarray as ar,exp,constants
from scipy.special import erfc
import impact_lib
reload(impact_lib)
from impact_lib import *
import gc
from scipy.ndimage.filters import gaussian_filter


plt.style.use('pub.mplstyle')
fig_path = r'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs'
halfwidth = 3.5
fullwidth = 7.2
# qcolors = {'Q1':'C0', 'Q2':'C1', 'Q3':'C2', 'Q4':'C3'}
# qcolors = {'Q1':'C0', 'Q2':'C4', 'Q3':'C3', 'Q4':'C1'}
# qcolors = {'Q1':'firebrick', 'Q2':'lightcoral', 'Q3':'mediumblue', 'Q4':'cornflowerblue'}
qcolors = {'Q1':'blue', 'Q2':'deepskyblue', 'Q3':'red', 'Q4':'magenta'}
run_plots = [
            # 'charge_jump',
            # 'fit_dropout',
            # 'gamma',
            # 'all_triggered',
            'offset_large',
            # 'offset_zoom',
            # 'qubit_spec',
            # 'ramsey_fit',
            # 'tracks',
            # 'pdfs',
            # 'hist1d_jumps',
            # 'hist2d_jumps_meas',
            # 'hist2d_jumps_sim',
            # 'bloch',
            # 'induced_rotation',
            # 'err_bitflip',
            # 'hist2d_err_phase',
            # 'err_phase_joint',
            # 'dipole',
            # 'spectra_emitted',
            # 'spectra_absorbed'
            ]

def draw_thresh_lines(ax, thresh=0.1):
    line_x = [-0.5,-thresh,np.nan,thresh,0.5]
    line_y = np.full(5,thresh)
    ax.plot( line_x,  line_y, 'k:' )
    ax.plot( line_x, -line_y, 'k:' )
    ax.plot(-line_y,  line_x, 'k:' )
    ax.plot( line_y,  line_x, 'k:' )

def no_tick_labels(*args):
    """Removes all tick labels from a vector of axes."""
    axs = np.array(args)
    for ax in np.ravel(axs):
        ax.set_xticklabels([])
        ax.set_yticklabels([])

def lower_left_tick_labels(*args, **kwargs):
    """Removes tick labels on all plots in vector of plots except those on the
    plot in the lower left."""
    axs = np.squeeze(np.array(args))
    if len(axs.shape) == 1:
        ni, nj = 1, axs.size
        axs = np.array([axs])
    elif len(axs.shape) == 2:
        ni, nj = axs.shape
    for i in range(ni):
        for j in range(nj):
            ax = axs[i,j]
            if (i == (ni - 1) and j == 0):
                if 'xlabel' in kwargs:
                    ax.set_xlabel(kwargs['xlabel'])
                if 'ylabel' in kwargs:
                    ax.set_ylabel(kwargs['ylabel'])
            else:
                ax.set_xticklabels([])
                ax.set_yticklabels([])

def edges_tick_labels(*args, **kwargs):
    """Removes tick labels from all plots in vector of plots but keeps the
    labels on the bottom and left of the array."""
    axs = np.squeeze(np.array(args))
    if len(axs.shape) == 1:
        ni, nj = 1, axs.size
        axs = np.array([axs])
    elif len(axs.shape) == 2:
        ni, nj = axs.shape
    for i in range(ni):
        for j in range(nj):
            ax = axs[i,j]
            if not i == (ni - 1):
                ax.set_xticklabels([])
            elif 'xlabel' in kwargs:
                ax.set_xlabel(kwargs['xlabel'])
            if not j == 0:
                ax.set_yticklabels([])
            elif 'ylabel' in kwargs:
                ax.set_ylabel(kwargs['ylabel'])

def format_hist2d(i, ax, cprofile, log=False, range=None, title=True):

    # color axes
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['left'].set_linewidth(3)
    q1,q2 = [('Q3','Q4'),('Q1','Q2'),('Q1','Q3')][i]
    ax.spines['bottom'].set_color(cprofile[q1])
    ax.spines['left'].set_color(cprofile[q2])
    
    ax.set_aspect(1)
    if range is not None:
        ax.set_xlim(range[0], range[1])
        ax.set_ylim(range[0], range[1])
        
    if title:
        title_str = ( u'$340\ \mathrm{\mu m}$',
                      u'$640\ \mathrm{\mu m}$',
                     u'$3195\ \mathrm{\mu m}$')
        ax.set_title(title_str[i])
    
    if log == True:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([1e-15,1e-12,1e-9,1e-6,1e-3])
        ax.set_yticks([1e-15,1e-12,1e-9,1e-6,1e-3])
        if not range:
                ax.set_xlim(1e-15, 1e-3)
                ax.set_ylim(1e-15, 1e-3)
        # ax.xaxis.labelpad = 0
        # ax.yaxis.labelpad = -6
    else:
        if not range:
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(-0.5, 0.5)
        ax.set_xticks([-0.5,-0.25,0,0.25,0.5])
        ax.set_yticks([-0.5,-0.25,0,0.25,0.5])
        ax.set_xticklabels([u'\u22120.5','',0,'',0.5])
        ax.set_yticklabels([u'\u22120.5','',0,'',0.5])
        ax.xaxis.labelpad = 0
        ax.yaxis.labelpad = -12
        

with open('dump_T1_sums.dat', 'rb') as f:
    P10,P01,P00,P11,n_trace,n0,n1,M1_before_trig = pickle.load(f)

"""Plot single charge jump trace"""
def plot_jump_trace(ax):
    # 2166.7, 1069.3, 1830.0, 1256.5, *497.8*, 73.8, 2104.8
    # for k,f in enumerate(files[:]): # 290
        # path, num = noiselib.path_to_num(f)
        # if num == 497:
            # print f
    p = (r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q4Corr\General\08-04-20\Charge_resetting\MATLABData\Charge_resetting_497.mat')
    data = noiselib.loadmat(p)
    charge_trace = data['Single_Shot_Occupations_RO2_SB1']
    # ax.set_xlabel('Time From Trigger')
    ax.set_ylabel('$P_{1}$')
    trig = 1570
    t = 0.040 * (np.arange(8000)-trig)
    ax.plot(t, noiselib.movingmean(charge_trace[8], 30), color=qcolors['Q1'] )
    # ax.arrow(0,1,0,-0.2, width=0.1, head_width=0.300, head_length=0.05, fc='black',
            # length_includes_head=True)
    # ax.plot([0],[noiselib.movingmean(charge_trace[8], 30)[trig]], 'xr', markersize=5)
    ax.set_xlim(-5,5)
    ax.axvline(x=0, color='k', linestyle=':')
if 'charge_jump' in run_plots:
    fig, ax = plt.subplots(1,1,figsize=(halfwidth,2))
    plot_jump_trace(ax)
    fig.savefig(fig_path+'\charge_jump.pdf')
    
"""Plot T1 dropouts in occupation."""
def plot_fit_occ(ax):
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

    N = 10000
    t = 0.04*np.arange(-N,N+1)
    ax.set_xlim(-5,5)
    ax.set_ylim(0.3,0.799)
    # ax.set_xlabel('Time From Trigger (ms)')
    ax.set_ylabel('$P_{1}$')
    for i,q in enumerate(['Q2','Q4']):
        o = noiselib.movingmean(1.*P01[q]/n0[q], 1)
        print n0[q][10000]
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
        print pcov
        ax.plot( t, o, label='Averaged P01 smoothed 5 {}'.format(q), color=qcolors[q] )
        ax.axvline(x=0, color='k', linestyle=':')
        if q == 'Q4':
            ax.plot( t, gausexp2(t,*popt), 'k', label='Averaged P01 smoothed 5 {}'.format(q) )
        # ax.plot( np.concatenate([t,t]), gausexp3(np.concatenate([t,t]),*popt), label='Averaged P01 smoothed 5 {}'.format(q) )
if 'fit_dropout' in run_plots:
    fig, ax = plt.subplots(1,1,figsize=(halfwidth,2))
    plot_fit_occ(ax)
    fig.savefig(fig_path+'\T1_occupation.pdf')

"""Plot T1 dropouts in T1 time."""
def plot_gamma(ax):
    axr = ax.twinx()
    ax.set_xlabel('Time from trigger (ms)')
    # ax.set_ylabel('Inferred T1 Time ($\mu$s)')
    ax.set_ylabel(r'$\Delta\Gamma_{01}$ ($\mathrm{\mu s}^{-1}$)')
    axr.set_ylabel(r'$\Delta x_{\mathrm{QP}}$ ($\times 10^{-6}$)')
    hbar = constants.physical_constants['Planck constant over 2 pi in eV s'][0]
    gap = 480e-6
    N = 10000
    t = 0.04*np.arange(-N,N+1)
    # a = {'Q2':0.6551, 'Q4':0.708}
    # c = {'Q2':0.1785, 'Q4':0.06842}
    # a = {'Q2':0.5385, 'Q4':0.5166}
    # c = {'Q2':0.2194, 'Q4':0.232}
    # a = {'Q2':0.7592, 'Q4':0.5423} # average a bunch of T1 curves and fit
    # c = {'Q2':0.1439, 'Q4':0.1198}
    a = {'Q2':1-2*0.38, 'Q4':1-2*0.31} # amp
    c = {'Q2':0.38, 'Q4':0.31} # offset
    f01 = {'Q2':(4.436+4.4308)/2, 'Q4':(4.4053+4.3995)/2}
    for i,q in enumerate(['Q2','Q4']):
        o = noiselib.movingmean(1.*P01[q]/n0[q], 1)
        T1 = -10./np.log((o-c[q])/a[q])
        gamma = 1./T1
        gamma = gamma - np.nanmean(gamma)
        gamma_to_xqp = np.pi/np.sqrt( 2*gap/hbar*2*np.pi*f01[q] ) / 1e-6
        ax.plot( t, gamma, label='$P_{}$ {}'.format('01',q), color=qcolors[q] )
        # axr.plot(t, gamma*gamma_to_xqp)
    ax.set_xlim(-5,5)
    ax.set_ylim(-0.1,0.38)
    axr.set_ylim(-0.1*gamma_to_xqp,0.38*gamma_to_xqp)
    ax.axvline(x=0, color='k', linestyle=':')
if 'gamma' in run_plots:
    fig, ax = plt.subplots(1,1,figsize=(halfwidth,2))
    plot_gamma(ax)
    fig.savefig(fig_path+'\T1_T1.pdf')

""" Plot all the above together """
if 'all_triggered' in run_plots:
    fig, (ax1,ax2,ax3) = plt.subplots(3, 1, figsize=(halfwidth,4.5))
    plot_jump_trace(ax1)
    plot_fit_occ(ax2)
    plot_gamma(ax3)
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    fig.savefig(fig_path+r'\all_triggered_traces.pdf')

"""Plot Charge Offset Large"""
def plot_charge_offset_large(ax):
    charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
    CO = ChargeOffset()
    for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
        CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
    # CO.plot_charge_offset()
    time, offset = CO.get_charge_offset()
    t0 = 200
    time = (time - time[t0])/60./60.
    for q in ['Q1','Q2','Q3','Q4']:
        offset[q] = offset[q] - offset[q][t0]
        ax.plot( time, offset[q], label=q, color=qcolors[q])#, linewidth=1. )
    ax.set_xlim(time[t0], time[t0]+10)
    trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
    y_max = np.vstack(offset.values())[:,trange].max()
    y_min = np.vstack(offset.values())[:,trange].min()
    ax.set_ylim(y_min, y_max)
    rect = patches.Rectangle((1,1.5),0.75,1.5,linewidth=1,edgecolor='k',facecolor='none')
    ax.add_patch(rect)
    rect = patches.Rectangle((175./60,-1),20./60,2,linewidth=1,edgecolor='k',facecolor='none')
    ax.add_patch(rect)
    # ax.legend()
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Offset Charge ($e$)')
if 'offset_large' in run_plots:
    fig, ax = plt.subplots(1, 1, figsize=(5.25,2.5))
    plot_charge_offset_large(ax)
    fig.savefig(fig_path+'\offset_charge.pdf')

"""Plot Charge Offset Zoom"""
def plot_charge_offset_zoom(ax1,ax2):
    charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
    CO = ChargeOffset()
    for d in [ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5]:
        CO.add_dataset(charge_path.format('fluxNoise', d['Q'], d['path_charge']))
    # CO.plot_charge_offset()
    time, offset = CO.get_charge_offset()
    t0 = 200
    time = (time - time[t0])/60.
    for q in ['Q1','Q2','Q3','Q4']:
        offset[q] = offset[q] - offset[q][t0]
        ax1.plot( time, offset[q], label=q, color=qcolors[q])#, linewidth=1. )
        ax2.plot( time, offset[q], label=q, color=qcolors[q])#, linewidth=1. )
    ax1.set_xlim(60, 105)
    ax2.set_xlim(175, 195)
    trange = np.where( (time>time[t0]) & (time<(time[t0]+10)) )[0]
    y_max = np.vstack(offset.values())[:,trange].max()
    y_min = np.vstack(offset.values())[:,trange].min()
    ax1.set_ylim(1.5, 3)
    ax2.set_ylim(-1, 1)
    # ax.legend()
    ax2.set_xlabel('Time (minutes)')
    # ax.set_ylabel('Charge Offset (e)')
if 'offset_zoom' in run_plots:
    fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(1.95,2.5))
    plot_charge_offset_zoom(ax1,ax2)
    fig.savefig(fig_path+'\offset_charge_zoom.pdf')

"""Plot Qubit Spec"""
def plot_qubit_spec(ax):
    path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1\General\07-19-20\Qubit_spectroscopy\MATLABData\Qubit_spectroscopy_001.mat'
    data = noiselib.loadmat(path)
    offset, vis, noff, nper = 0, 4.1, .2, 1.05*0.2709
    f_qb = (data['Qubit_Frequency']  - 4.56415)*1000
    df_qb = f_qb[1] - f_qb[0]
    bias = data['Qubit_Slow_Bias_1']/nper
    dbias = bias[1] - bias[0]
    o = np.fliplr(data['Single_Shot_Occupation_SB1'].T)
    p = ax.imshow(o, origin='lower', aspect='auto', cmap='gray_r',
                  extent=( bias[0]-dbias/2., bias[-1]+dbias/2.,
                           f_qb[0]-df_qb/2., f_qb[-1]+df_qb/2.) )
    o = offset + vis/2*np.cos((np.pi*bias*1.05-noff)/1.05)
    ax.plot(bias,  o, 'C0', linewidth=2)
    ax.plot(bias, -o, 'C3', linewidth=2)
    ax.set_xlabel('Applied offset charge ($e$)')
    ax.set_ylabel('$f-\overline{f_{01}}$ (MHz)')
if 'qubit_spec' in run_plots:
    fig, ax = plt.subplots(1, 1, figsize=(2.55,2))
    plot_qubit_spec(ax)
    fig.savefig(fig_path+'\qubit_spec.pdf')

"""Plot Ramsey Fit"""
def plot_ramsey_fit(ax):
    path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Q4Corr\General\08-30-20\Charge_resetting\MATLABData\\'
    # 23-24, 28-29, 34-35, 64-65, 71-73
    data1 = noiselib.loadmat(path+'Charge_resetting_024.mat')
    data2 = noiselib.loadmat(path+'Charge_resetting_023.mat')
    o1 = np.flip(data1['Single_Shot_Occupation_RO2_SB1'])
    o2 = np.flip(data2['Single_Shot_Occupation_RO2_SB1'])
    b = np.arange(-0.5,0.5,0.01)
    break_index = 48 - 1
    guess_offset = np.mean(o1)
    guess_vis = np.max(o1) - np.min(o1)
    def ram(bias, offset, vis, noff, nper):
        return offset + vis/2.*np.cos(np.pi*np.cos(np.pi*(bias-noff)/nper))
    popt1, pcov1 = curve_fit(ram, b, o1, 
                                p0=[guess_offset,guess_vis,-0.05,0.3])
    popt2, pcov2 = curve_fit(ram, b[b<b[break_index]], o1[b<b[break_index]], 
                                p0=[guess_offset,guess_vis,0.16,0.3])
    print popt1,popt2
    popt1[0] = popt1[0] - 0.02
    popt2p = popt1.copy()
    popt2p[2] = 0.1529 #popt2[2]
    ax.plot(b/popt1[-1], o1, 'C1.', markersize=6)
    ax.plot(b/popt1[-1], o2, 'C2.', markersize=6)
    ax.plot(b/popt1[-1], ram(b,*popt1), 'C1', linewidth=1.5)
    ax.plot(b[break_index:]/popt1[-1], ram(b[break_index:], *popt2p), 'C2', linewidth=1.5)
    ax.plot(b[:break_index+1]/popt1[-1], ram(b[:break_index+1], *popt2p), 'C2:', linewidth=1.5)
    ax.set_xlabel('Applied offset charge ($e$)')
    ax.set_ylabel('$P_1$')
    ax.set_yticks([0.2,0.4,0.6,0.8])
if 'ramsey_fit' in run_plots:
    fig, ax = plt.subplots(1, 1, figsize=(2.55,2))
    plot_ramsey_fit(ax)
    fig.savefig(fig_path+r'\ramsey_fit.pdf')

"""Plot Muon and Gamma tracks"""
def plot_tracks(ax1,ax2,ax3,ax4):
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
    _range = {'Muons':range(3,33),'Gamma':range(0,60)}
    for k,hit_type in enumerate(['Muons','Gamma']):
        app = Controller(sys.argv, fQ=0.01, calc_all=False, plot=False,
                            pdfs_file='{}/ChargePDFs_{}r_vhs_fine.npy'.format(path,500),
                            event_files=[path+'/{}.txt'.format(hit_type)])
        track_xyz, track_color, charge_xyz, charge_polarity = [],[],[],[]
        for i in _range[hit_type]:
            app.i = i
            track_xyz_, track_color_, charge_xyz_, charge_polarity_ = app.sim_event()
            track_xyz       += [np.full((1,3),np.nan)] + [track_xyz_]
            track_color     += [np.full((1,3),np.nan)] + [track_color_]
            charge_xyz      += [np.full((1,3),np.nan)] + [charge_xyz_]
            charge_polarity += [np.full((1,3),np.nan)] + [charge_polarity_]
        track_xyz = np.concatenate(track_xyz)
        # charge_xyz = np.concatenate(charge_xyz)
        # charge_polarity = np.concatenate(charge_polarity)
        del app
        gc.collect()
        ax1.scatter(track_xyz[:,0], track_xyz[:,1], s=1, c=['C0','C3'][k])
        ax3.scatter(track_xyz[:,0], track_xyz[:,2], s=1, c=['C0','C3'][k])
        ax2.scatter(track_xyz[:,2], track_xyz[:,1], s=1, c=['C0','C3'][k])
        if hit_type == 'Muons':
            n = 7
            for ax,i,j in [(ax1,0,1),(ax3,0,2),(ax2,2,1)]:
                pos = charge_xyz[n][charge_polarity[n]>0]
                neg = charge_xyz[n][charge_polarity[n]<0]
                alli = np.vstack([pos[:,i],neg[:,i]]).reshape(-1,order='F')
                allj = np.vstack([pos[:,j],neg[:,j]]).reshape(-1,order='F')
                ax.scatter(alli, allj, s=0.25, c=('C1','C2'), zorder=5)
        qubit_xys = [( 1.564,  0.570),
                     ( 1.564, -0.070),
                     (-1.564, -0.080),
                     (-1.564, -0.420)]
        for x,y in qubit_xys:
            ax1.add_artist( plt.Circle((x,y), 0.07, color='black', zorder=10) )
            ax1.add_artist( plt.Circle((x,y), 0.095, ec='black', zorder=10, 
                                        fill=False, linewidth=0.5) )
    # ax1.set_aspect(1)
    # ax2.set_aspect(1./2.5)
    # ax3.set_aspect(2.5)
    ax1.set_xlim(-6.25/2,6.25/2)
    ax1.set_ylim(-6.25/2,6.25/2)
    ax2.set_ylim(-6.25/2,6.25/2)
    ax2.set_xlim(-0.375/2,0.375/2)
    ax3.set_xlim(-6.25/2,6.25/2)
    ax3.set_ylim(-0.375/2,0.375/2)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax3.set_xticks([])
    ax4.set_yticks([])
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
if 'tracks' in run_plots:
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2, figsize=(halfwidth,halfwidth), 
                                            gridspec_kw={'width_ratios': [1, 7],
                                                         'height_ratios': [1, 7]})
    plot_tracks(ax4,ax3,ax2,ax1)
    fig.delaxes(ax1)
    fig.savefig(fig_path+'\hit_tracks.pdf')

"""Plot PDFs"""
def plot_pdf_raw(ax1,ax2):
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
    PDFs = np.load(path+'/ChargePDFs_300r_vhs_fine.npy',allow_pickle=True).tolist()
    e = ImpactEvent([PDFs], np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
    
    for j,species in enumerate(['electrons','holes']):
        
        zcuts = e.PDFs[species]['zarray']
        zind = np.argmin(np.abs(0-zcuts))
        xyz = e.PDFs[species][zcuts[zind]][species]
        
        xyz[51:,:,:] += np.flip(xyz[:50,:,:], axis=0)
        xyz[51:,51:,:] += np.flip(xyz[51:,:50,:], axis=1)
        xyz[51:,51:,51:] += np.flip(xyz[51:,51:,:50], axis=2)
        xyz[51:,51:,51:] /= 8.
        
        z = e.PDFs[species]['z']
        inds = (z < 0.025) & (z > -0.0000025)
        xy = np.sum(xyz[:,:,inds],axis=-1)
        
        ax = [ax1,ax2][j]
        ax.imshow(xy, origin='lower', 
                  norm=mpl.colors.LogNorm(), 
                  cmap=cm.get_cmap('jet'))

        # Y, X = np.meshgrid((x[1:]+x[0:-1])/2,(y[1:]+y[0:-1])/2)
        # ax.contour(X, Y, np.log10(h), levels=[2.5,3,3.5,4,4.5,5,5.5,6], color='black')

        # ax3.contour(gaussian_filter(np.sum(v00.T, axis=2),0.9), norm=mpl.colors.LogNorm(), origin='lower')
            # # ax1.contour(v00.T[:,:,25], norm=mpl.colors.LogNorm(), origin='lower')
    
        # ax.set_xlim(-0.1,0.1)
        # ax.set_ylim(-0.1,0.1)
        # ax.set_xticks([-.1,-0.05,0,0.05,.1])
        # ax.set_xticklabels([-.1,'',0,'',.1])
        ax.set_xticks([])
        ax.set_xticklabels([])
        # ax.set_yticks([-.1,-0.05,0,0.05,.1])
        # ax.set_yticklabels([-.1,'',0,'',.1])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
def plot_pdf_cut(ax1,ax2):
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
    PDFs = np.load(path+'/ChargePDFs_300r_vhs_fine.npy',allow_pickle=True).tolist()
    e = ImpactEvent([PDFs], np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
    
    r110, r100 = {}, {}
    n110, n100 = {}, {}
    
    for j,species in enumerate(['electrons','holes']):
        
        zcuts = e.PDFs[species]['zarray']
        zind = np.argmin(np.abs(0-zcuts))
        xyz = e.PDFs[species][zcuts[zind]][species]
        
        # this code flips the quadrants into the +++ quadrant and averages it
        # xyz[51:,:,:] += np.flip(xyz[:50,:,:], axis=0)
        # xyz[51:,51:,:] += np.flip(xyz[51:,:50,:], axis=1)
        # xyz[51:,51:,51:] += np.flip(xyz[51:,51:,:50], axis=2)
        # xyz[51:,51:,51:] /= 8.
        
        x,y,z = e.PDFs[species]['x'], e.PDFs[species]['y'], e.PDFs[species]['z']
        r110[species] = e.PDFs[species]['x']
        r100[species] = e.PDFs[species]['x'] * np.sqrt(2)
        n110[species] = xyz[:,50,50]
        n100[species] = xyz[np.arange(101),np.arange(101),50]
        # normalize
        # n110[species] = n110[species] / np.sum(n110[species])
        # n100[species] = n100[species] / np.sum(n100[species])
        
        ax = [ax1,ax2][j]
        ax.plot(r110[species], n110[species])
        ax.plot(r100[species], n100[species])
        ax.plot(r110[species], n110[species].max() * np.exp(-r110[species]/300e-3), 'k')
        ax.set_yscale('log')
        
        # average all points with the same radius
        xi, yi, zi = np.indices(xyz.shape) - np.full(xyz.shape, 50)
        dx,dy,dz = e.PDFs[species]['dx'], e.PDFs[species]['dy'], e.PDFs[species]['dz']
        ri = np.sqrt((xi*dx)**2 + (yi*dy)**2 + (zi*dz)**2)
        # ri = np.sqrt(x**2 + y**2 + z**2)
        i_sorted = np.argsort(np.ravel(ri))
        r_sorted = np.ravel(ri)[i_sorted]
        n_sorted = np.ravel(xyz)[i_sorted]
        r_avg, n_avg = [], []
        r0 = r_sorted[0]
        n0 = []
        for i in range(r_sorted.size):
            if r_sorted[i] == r0:
                n0 += [n_sorted[i]]
            else:
                r_avg += [r0]
                n_avg += [np.mean(n0)]
                r0 = r_sorted[i]
                n0 = [n_sorted[i]]
        # ax.plot(r_sorted, n_sorted)
        ax.plot(r_avg, n_avg)
        
    fig, ax = plt.subplots(1, 1)
    ax.plot(r110['electrons'], n110['holes']-n110['electrons'])
    ax.plot(r100['electrons'], n100['holes']-n100['electrons'])
def plot_pdf_qgen(ax1,ax2):
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
    PDFs = np.load(path+'/ChargePDFs_300r_vhs_fine.npy',allow_pickle=True).tolist()
    e = ImpactEvent([PDFs], np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
    res = e.diffuseCharge(0.0, 0.0, 0.0, 1e8)#, PDFdict['fine'])
    
    # # plot for (electrons - holes)
    # fig, axxx = plt.subplots(1, 1, figsize=(halfwidth,halfwidth))
    # inds = (res['electrons'][2] < 0.025) & (res['holes'][2] > -0.025)
    # h,x,y,i = axxx.hist2d(res['electrons'][0][inds] - res['holes'][0][inds],
                         # res['electrons'][1][inds] - res['holes'][1][inds],
                         # bins=(101,101), range=[[-0.25,0.25],[-0.25,0.25]],
                         # norm=mpl.colors.LogNorm(vmin=1e0,vmax=1e6),
                         # cmap=cm.get_cmap('jet'))
    # Y, X = np.meshgrid((x[1:]+x[0:-1])/2,(y[1:]+y[0:-1])/2)
    # axxx.contour(X, Y, np.log10(h), levels=[2.5,3,3.5,4,4.5,5,5.5,6], color='black')

    for j,species in enumerate(['electrons','holes']):
    
        ax = [ax1,ax2][j]

        inds = (res[species][2] < 0.025) & (res[species][2] > -0.025)

        h,x,y,i = ax.hist2d(res[species][0][inds],
                         res[species][1][inds],
                         bins=(101,101), range=[[-0.25,0.25],[-0.25,0.25]],
                         norm=mpl.colors.LogNorm(vmin=1e0,vmax=1e6),
                         cmap=cm.get_cmap('jet'))

        Y, X = np.meshgrid((x[1:]+x[0:-1])/2,(y[1:]+y[0:-1])/2)
        ax.contour(X, Y, np.log10(h), levels=[2.5,3,3.5,4,4.5,5,5.5,6], color='black')

        ax.set_xlim(-0.1,0.1)
        ax.set_ylim(-0.1,0.1)
        # ax.set_xticks([-.1,-0.05,0,0.05,.1])
        # ax.set_xticklabels([-.1,'',0,'',.1])
        ax.set_xticks([])
        ax.set_xticklabels([])
        # ax.set_yticks([-.1,-0.05,0,0.05,.1])
        # ax.set_yticklabels([-.1,'',0,'',.1])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
if 'pdfs' in run_plots:
    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(halfwidth,halfwidth/2.))
    # plot_pdf_raw(ax1,ax2)
    plot_pdf_qgen(ax1,ax2)
    # plot_pdf_cut(ax1,ax2)
    # fig.savefig(fig_path+'\pdfs.pdf')

""" Plot distribution of charge rotations induced on qubit. """
def plot_rotation_hist(ax, rot_induced, i):
    err = rot_induced
    h,bins = np.histogram(err, bins=100000)
    center = (bins[1:] + bins[:-1])/2.
    ax.step(center, 1.*h/np.sum(h), where='mid', zorder=-i, color=['C3','C0'][i])
    ax.set_yscale('log')
    ax.set_ylabel('Normalized\nCounts')
    ax.set_xlabel('Rotation on Bloch sphere (rad)')
    ax.set_xscale('log')
    ax.set_xlim(1e-4,1e-1)
    # ax.set_title('Distribution of induced rotation')
    plt.draw()
    plt.pause(0.05)
def plot_error_hist(ax, rot_induced):
    # err = rot_induced[500,0.1]**2/4.
    err = np.sin(rot_induced/2.)**2
    h,bins = np.histogram(err, bins=1500)
    center = (bins[1:] + bins[:-1])/2.
    ax.bar(center, 1.*h/np.sum(h), width=center[1]-center[0])
    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('Counts')
    ax.set_xlabel('Qubit Error')
    # ax.set_title('Distribution of induced rotation')
    plt.draw()
    plt.pause(0.05)
def plot_error(ax, rot_induced, cindex):
    # err = rot_induced[500,0.1]**2/4.
    err = np.sin(rot_induced/2.)**2
    # for i in range(4):
        # ax.plot(np.sort(err[:,i]), 1.*np.arange(err[:,i].size)/err[:,i].size,
                # linewidth=0.5, alpha=0.5, color=qcolors['Q{}'.format(i+1)])
    sorted_err = np.sort(np.ravel(err))
    ax.plot(sorted_err, 1.*np.arange(sorted_err.size)/sorted_err.size, #'.')
            linewidth=2, color=['C3','C0'][i])
    ax.axvline(x=1e-3, color='k', linestyle=':')
    ax.axvline(x=1e-4, color='k', linestyle=':')
    y3 = 1.*np.searchsorted(sorted_err, 1e-3)/sorted_err.size
    y4 = 1.*np.searchsorted(sorted_err, 1e-4)/sorted_err.size
    print 1.-y3,1.-y4
    ax.plot([0,1e-3],[y3,y3],':', color=['C3','C0'][i])
    ax.plot([0,1e-4],[y4,y4],':', color=['C3','C0'][i])
    ax.set_xlim(1e-10,1)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylabel('Normalized\nIntegrated Counts')
    ax.set_xlabel('Qubit Error')
    # ax.set_title('Distribution of induced rotation')
    plt.draw()
    plt.pause(0.05)
if 'induced_rotation' in run_plots:
    fig, (ax_rot,ax_err) = plt.subplots(2, 1, figsize=(halfwidth,4))
    for i,hit_type in enumerate(['gammas','muons']):
        with open('dump_sim_impacts_{}_vhs.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        rot_induced = rot_induced[300,0.2]
        plot_rotation_hist(ax_rot, rot_induced, i)
        plot_error(ax_err, rot_induced, i)
    fig.savefig(fig_path+r'\err.pdf')
    
def plot_hist1d_jumps_meas(ax, q):
    h, bins = np.histogram(q, bins=500, range=(-0.5,0.5))
    x = (bins[1:]+bins[:-1])/2
    center = (bins[:-1] + bins[1:])/2
    center2 = center**2 * np.sign(center)
    bg = 1. * np.sum(np.abs(q)>0.05) / np.sum(h)
    # ax.step(center, 6333./1./bg*noiselib.movingmean(h,3), where='mid' )
    ax.step(center, 1.*noiselib.movingmean(h,3)/np.sum(h), where='mid' )
    # ax.set_ylim([10e-1, 1.5*h.max()])
    return bg
def plot_hist1d_jumps_sim(ax, q_induced, bg_meas):
    q = q_induced[300,0.2][:,1-1]
    r = np.random.normal(0, 0.0132688, q.size)
    e = noiselib.alias(q+r)
    h, bins = np.histogram(e, bins=500, range=(-0.5,0.5))
    x = (bins[1:]+bins[:-1])/2
    center = (bins[:-1] + bins[1:])/2
    center2 = center**2 * np.sign(center)
    bg_sim = 1. * np.sum(np.abs(e)>0.05) / np.sum(h)
    # ax.step(center, 6333./1./bg*noiselib.movingmean(h,3), where='mid' )
    ax.step(center, 1.*noiselib.movingmean(h,3)/np.sum(h)*bg_meas/bg_sim, where='mid' )
def set_style(ax):
    ax.set_title('')
    ax.set_xlabel('$\Delta q$ ($e$)')
    ax.set_ylabel('Normalized Counts')
    ax.set_xlim(-0.5,0.5)
    ax.set_yscale('log')
if 'hist1d_jumps' in run_plots:
    fig, ax = plt.subplots(1, 1, figsize=(fullwidth,3))
    # fig = plt.figure(700)
    # ax = fig.axes[0]
    with open('dump_measured_Q1.dat'.format(hit_type), 'rb') as f:
        q = pickle.load(f)
        bg = plot_hist1d_jumps_meas(ax, q)
    for i,hit_type in enumerate(['gammas']):
        with open('dump_sim_impacts_{}_vhs.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        plot_hist1d_jumps_sim(ax, q_induced, bg)
    set_style(ax)
    fig.savefig(fig_path+r'\hist1d_jumps.pdf')
    
    
def plot_hist2d(axs, q, title=True, log=False):
    ij = [(2,3),(0,1),(0,2)]
    q = noiselib.alias(q, 0.5)
    for j,ax in enumerate(axs):
        ax.scatter(q[:,ij[j][0]], q[:,ij[j][1]], c='k', marker='.', s=4)
        format_hist2d(j, ax, cprofile=qcolors, title=title, log=log)
def add_noise(q):
    sigma = [[0.02086807133980249, 0.04472867181129812, 
              0.008811574164510916, 0.011038152654827339]]
    m_sigma = np.repeat(sigma, q.shape[0], axis=0)
    return q + np.random.normal(0, m_sigma)
if 'hist2d_jumps_meas' in run_plots:
    fig, axs = plt.subplots(1, 3, figsize=(fullwidth,2.6))
    with open('dump_measured_Q1234.dat'.format(hit_type), 'rb') as f:
        q = pickle.load(f)
    plot_hist2d(axs, q)
    for ax in axs:
        draw_thresh_lines(ax)
    lower_left_tick_labels(axs, xlabel='$\Delta q_\mathrm{1}$ ($e$)',
                                ylabel='$\Delta q_\mathrm{2}$ ($e$)')
    plt.draw()
    plt.pause(0.05)
    fig.savefig(fig_path+r'\hist2d_jumps_meas.pdf')
if 'hist2d_jumps_sim' in run_plots:
    fig, axs = plt.subplots(2, 3, figsize=(fullwidth-0.5,5.2-0.5))
    for i,hit_type in enumerate(['gammas','muons']):
        with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        # q_induced = add_noise(q_induced)
        plot_hist2d(axs[i,:], q_induced[300,0.2], title=(i==0))
    lower_left_tick_labels(axs, xlabel='$\Delta q_\mathrm{1}$ ($e$)',
                                ylabel='$\Delta q_\mathrm{2}$ ($e$)')
    fig.savefig(fig_path+r'\hist2d_jumps_sim.pdf')
if 'hist2d_err_phase' in run_plots:
    fig, axs = plt.subplots(2, 3, figsize=(fullwidth-0.5,5.2-0.5))
    for i,hit_type in enumerate(['gammas','muons']):
        with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        dq = q_induced[300,0.2]
        err = dw01**2 / 12. * np.sin(np.pi*dq/2.)**2 * tau**2
        err = err[~np.all(err == 0, axis=1)]
        plot_hist2d(axs[i,:], err, title=(i==0), log=True)
    edges_tick_labels(axs, xlabel=u'$\epsilon_{\phi,1}$', 
                           ylabel=u'$\epsilon_{\phi,2}$')
    fig.savefig(fig_path+r'\hist2d_err_phase.pdf')
     
def plot_err_phase_joint(ax, q_induced):
    ij = [(2,3),(0,1),(0,2)]
    for j in range(3):
        err = dw01**2 / 12. * np.sin(np.pi*q_induced/2.)**2 * tau**2
        err = err[~np.all(err == 0, axis=1)]
        err_axis = np.logspace(-15, -3, 101)
        joint_err = [ np.sum( (err[:,ij[j][0]]>t)&(err[:,ij[j][1]]>t) )
                      for t in err_axis ]
        joint_err = 1. - 1. * np.array(joint_err) / err.shape[0]
        ax.plot( err_axis, joint_err )
    ax.set_xscale('log')
    ax.set_xlim(1e-15,1e-3)
    ax.set_xticks([1e-15,1e-12,1e-9,1e-6,1e-3])
if 'err_phase_joint' in run_plots:
    fig, axs = plt.subplots(1, 2, figsize=(fullwidth,3))
    for i,hit_type in enumerate(['gammas','muons']):
        with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        plot_err_phase_joint(axs[i], q_induced[300,0.2])
    edges_tick_labels(axs, xlabel=u'$\epsilon_\phi$', ylabel='Normalized Integrated Counts')
    fig.savefig(fig_path+r'\err_phase_joint.pdf')


def map_alpha(r, z, q, rp, zp): 
    # my own bilinear interpolation, much faster
    # https://math.stackexchange.com/questions/3230376/
    #   interpolate-between-4-points-on-a-2d-plane
    ri = np.searchsorted(r, rp)
    zi = np.searchsorted(z, zp)
    dr_, dz_ = r[1]-r[0], z[1]-z[0]
    rpp = (rp-r[np.clip(ri-1,0,r.size-1)])/dr_
    zpp = (zp-z[np.clip(zi-1,0,z.size-1)])/dz_
    q1 = q[np.clip(zi-1,0,z.size-1), np.clip(ri-1,0,r.size-1)]
    q2 = q[np.clip(zi-1,0,z.size-1), np.clip(ri  ,0,r.size-1)]
    q3 = q[np.clip(zi  ,0,z.size-1), np.clip(ri  ,0,r.size-1)]
    q4 = q[np.clip(zi  ,0,z.size-1), np.clip(ri-1,0,r.size-1)]
    return (1-rpp)*(1-zpp)*q1 + rpp*(1-zpp)*q2 + (1-rpp)*zpp*q3 + rpp*zpp*q4
if 'err_bitflip' in run_plots:
    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(fullwidth,1.5))
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data/charge_map.mat'
    q_map = noiselib.loadmat(path)
    q = -q_map['charge_mat']
    q[~np.isfinite(q)] = -1. # right below center of qubit
    r, z = q_map['r']*1e-3, q_map['z']*1e-3
    dr, dz = np.diff(r)[0], np.diff(z)[0]
    dqr, dqz = np.gradient(q, dz*1e-3, dr*1e-3) # *1e-3 factor for units m instead of mm
    dq = np.sqrt(dqr**2 + dqz**2)
    p = ax1.imshow(q, extent=(r[0]-dr/2., r[-1]+dr/2., z[0]-dz/2., z[-1]+dz/2.),
                   norm=mpl.colors.LogNorm() )
    c = fig.colorbar(p, ax=ax1)
    p = ax2.imshow(dq, extent=(r[0]-dr/2., r[-1]+dr/2., z[0]-dz/2., z[-1]+dz/2.),
                   norm=mpl.colors.LogNorm() )
    c = fig.colorbar(p, ax=ax2)
    ax1.set_xlim(0,1)
    ax2.set_xlim(0,1)
    dalpha_map  = lambda rp,zp: map_alpha( r, z, dq,  rp, zp )
    dalphar_map = lambda rp,zp: map_alpha( r, z, dqr, rp, zp )
    dalphaz_map = lambda rp,zp: map_alpha( r, z, dqz, rp, zp )
    
    qubit_xys = [( 1.564,  0.570),
                 ( 1.564, -0.070),
                 (-1.564, -0.080),
                 (-1.564, -0.420)]
                 
    cs = 1e4
    # w0 = 2 * np.pi * 4.1e9
    # Ec = 2 * np.pi * 435e6
    w0 = 2 * np.pi * 5e9
    Ec = 2 * np.pi * 300e6
    
    fig1d, axs1d = plt.subplots(1, 2, figsize=(fullwidth,2.5))
    fig2d, axs2d = plt.subplots(2, 3, figsize=(fullwidth-0.5,4.6))
    
    for h,hit_type in enumerate(['gammas','muons']):
        data = []
        path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
        if hit_type == 'gammas':
            efiles = [path+'/Gamma.txt', 
                      path+'/Gamma_10deg.txt',
                      path+'/Gamma_10deg_pt2.txt']
        elif hit_type == 'muons':
            efiles = [path+'/Muons.txt', path+'/Muons2.txt']
        for path in efiles:
            data.append( np.loadtxt(path, skiprows=1) )
        data = np.concatenate(data)
        split_data = np.split(data, np.where(np.diff(data[:,0]))[0]+1)
        
        E_sfq = []
        for event_data in split_data:
            xyz, E = event_data[:,[2,3,4]], event_data[:,5] * 1e6
            n = 2*np.floor(E/3.6)
            # pick random dipole angle
            phi = 2 * np.pi * np.random.random(E.size)
            costheta = np.random.uniform(-1, 1, E.size)
            theta = np.arccos(costheta)
            dip_x = np.sin( theta ) * np.cos( phi )
            dip_y = np.sin( theta ) * np.sin( phi )
            dip_z = np.cos( theta )
            E_sfq_q = []
            for i,xy_q in enumerate(qubit_xys):
                dx,dy,dz = (xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
                dr = np.hypot(dx,dy)
                da   = dalpha_map(dr,-dz)
                da_r = dalphar_map(dr,-dz)
                da_z = dalphaz_map(dr,-dz)
                da_x, da_y = np.cos(da_r/da), np.sin(da_r/da)
                E_sfq_q += [np.sum( n * cs**2 / w0**2 * ((da_x*dip_x)**2 
                                                       + (da_y*dip_y)**2 
                                                       + (da_z*dip_z)**2) )]
            E_sfq += [E_sfq_q]
        E_sfq = np.array(E_sfq)
        err = E_sfq * 2./3. * Ec/w0
        print('% above error of {{1e-4,1e-5,1e-6}} = {},{},{}'.format(
                    100.*np.sum(err>1e-4)/err.size, 
                    100.*np.sum(err>1e-5)/err.size,
                    100.*np.sum(err>1e-6)/err.size ))
        print( 100.*np.sum((err[:,3-1]>1e-8) & (err[:,4-1]>1e-8))/err.shape[0],
               100.*np.sum((err[:,1-1]>1e-8) & (err[:,2-1]>1e-8))/err.shape[0] )
        
        """
        th = np.linspace(0,np.pi,100)
        dist = 0.5 * np.sin(th)
        int_dist = np.cumsum(dist)/np.sum(dist)
        for event_data in split_data:
            xyz, E = event_data[:,[2,3,4]], event_data[:,5] * 1e6
            E_sfq_q = []
            for i,xy_q in enumerate(qubit_xys):
                dx,dy,dz = (xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
                dr = np.hypot(dx,dy)
                da = dalpha_map(dr,-dz)
                n = 2*np.floor(E/3.6)
                # pick a random angle
                rnd = np.random.random(E.size)
                inds = np.searchsorted(int_dist, rnd)
                theta = th[inds]
                E_sfq_q += [np.sum( n * cs**2 * da**2 / w0**2 * np.cos(theta)**2 )]
            E_sfq += [E_sfq_q]
        E_sfq = np.array(E_sfq)
        """
        
        hist, bins = np.histogram(np.mean(E_sfq,axis=1),
                      bins=np.logspace(np.log10(1e-14),np.log10(1e0),101) )
        axs1d[h].bar( bins[:-1], 1.*hist/np.sum(hist), width=np.diff(bins), align='edge' )
        axs1d[h].set_xlabel('$E/E_c$')
        axs1d[h].set_ylabel('Normalized Counts')
        axs1d[h].set_xscale('log')
        axs1d[h].set_yscale('log')
        # axs1d[h].set_title(hit_type)
        
        for i,ax in enumerate(axs2d[h]):
            j,k = [(3,4),(1,2),(1,3)][i]
            ax.scatter(err[:,j-1], err[:,k-1], c='k', marker='.', s=4)
            format_hist2d(i, ax, qcolors, log=True, range=(1e-15,1e-1), title=(h==0))
        
    edges_tick_labels(axs2d, xlabel='$\epsilon_1$', ylabel='$\epsilon_2$')
    
    fig1d.savefig(fig_path+r'\hist1d_err_bitflip.pdf')
    fig2d.savefig(fig_path+r'\hist2d_err_bitflip.pdf')


if 'spectra_absorbed' in run_plots:
    path = r'Z:\mcdermott-group\data\fluxNoise2\sim_data'
    E_muons, rate_muons = np.loadtxt(path+'\simulation_muons.txt', skiprows=1).T
    E_gammas, rate_gammas = np.loadtxt(path+'\simulation_ambient_radioactivity.txt', skiprows=1).T
    # change units to counts/sec/keV
    rate_muons /= (0.03*1000)
    rate_gammas /= (0.03*1000)
    
    fig, ax = plt.subplots(1, 1, figsize=(halfwidth,2.5))
    ax.step(E_gammas, rate_gammas, 'C3', where='mid')
    ax.step(E_muons, rate_muons, 'C0', where='mid')
    ax.set_xlim(0, 3)
    ax.set_yscale('log')
    ax.set_xlabel('Energy (MeV)')
    ax.set_ylabel('Rate (Counts/sec/keV)')
    fig.savefig(fig_path+r'\spectra_absorbed.pdf')
    
def get_calibrated_spect_lngs():
    return np.linspace(0,3000,1001), np.full(1001, 1e-4)
def get_calibrated_spect_madison():
    return np.linspace(0,3000,1001), np.full(1001, 1e-4)
if 'spectra_emitted' in run_plots:
    path = r'Z:\mcdermott-group\data\fluxNoise2\sim_data'
    E_mad, rate_mad = get_calibrated_spect_madison()
    E_lngs, rate_lngs = get_calibrated_spect_lngs()
    elem_labels = {'$^{214}\mathrm{Bi}$':600, '$^{40}\mathrm{K}$':1500, 
      '$^{214}\mathrm{Bi}$':1700, '$^{214}\mathrm{Bi}$':2200, '$^{208}\mathrm{Tl}$':2600}
    
    fig, ax = plt.subplots(1, 1, figsize=(halfwidth,2.5))
    ax.step(E_mad, rate_mad, 'C0', where='mid')
    ax.step(E_lngs, rate_lngs, 'C3', where='mid')
    # add labels
    for label in elem_labels:
        v_mad  =  rate_mad[np.searchsorted(E_mad, elem_labels[label])]
        v_lngs = rate_lngs[np.searchsorted(E_lngs,elem_labels[label])]
        ax.text( elem_labels[label], 2*max(v_mad,v_lngs), label, 
                 horizontalalignment='center')
    ax.set_xlim(0, 3000)
    ax.set_yscale('log')
    ax.set_xlabel('Energy (MeV)')
    ax.set_ylabel('Rate (Counts/sec/keV)')
    fig.savefig(fig_path+r'\spectra_emitted.pdf')


def plot_dipole_sphere(ax):
    # Create a sphere
    r = 1
    phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
    x = r*np.sin(phi)*np.cos(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(phi)
    ax.plot_surface(
        x, y, z,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)
    
    n = 200
    phi = 2 * np.pi * np.random.random(n)
    costheta = np.random.uniform(-1, 1, n)
    theta = np.arccos(costheta)
    xx = np.sin( theta ) * np.cos( phi )
    yy = np.sin( theta ) * np.sin( phi )
    zz = np.cos( theta )
    
    ax.scatter(xx[:n/2], yy[:n/2], zz[:n/2], color="C1", s=20, alpha=1.)
    ax.scatter(xx[n/2:], yy[n/2:], zz[n/2:], color="C3", s=20, alpha=1.)

    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("equal")
    plt.tight_layout()
if 'dipole' in run_plots:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_dipole_sphere(ax)
    

"""Plot Bloch Sphere"""
def circle(x=0, y=0, r=1., theta1=0, theta2=2*np.pi):
    theta = np.linspace(theta1,theta2,100)
    xp = x + r * np.cos(theta)
    yp = y + r * np.sin(theta)
    return xp,yp
if 'bloch' in run_plots:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_axis_off()
    xc,yc = circle()
    zc = np.full(100,0)
    ax.plot(xc, yc, zc, 'k--')
    ax.plot(xc, zc, yc, 'k')
    ax.plot(zc, xc, yc, 'k')

    # ax.quiver(-2, 0, 0, 1, 0, 0, color='k', length=4, headwidth=1)

plt.draw()
plt.pause(0.05)