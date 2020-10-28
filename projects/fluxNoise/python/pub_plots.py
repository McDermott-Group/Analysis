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
            # 'offset_large',
            # 'offset_zoom',
            # 'qubit_spec',
            # 'ramsey_fit',
            'tracks',
            # 'pdfs',
            # 'hist1d',
            # 'hist2d',
            # 'bloch',
            # 'induced_rotation',
            ]

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
    ax.arrow(0,1,0,-0.2, width=0.1, head_width=0.300, head_length=0.05, fc='black',
            length_includes_head=True)
    ax.plot([0],[noiselib.movingmean(charge_trace[8], 30)[trig]], 'xr', markersize=5)
    ax.set_xlim(-5,5)
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
    ax.set_ylim(0.3,0.8)
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
    ax.set_ylabel('Charge Offset ($e$)')
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
    fig, ax = plt.subplots(1, 1, figsize=(fullwidth/3.,2))
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
    fig, ax = plt.subplots(1, 1, figsize=(fullwidth/3.,2))
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
def plot_pdf(ax1,ax2):
    path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
    PDFs = np.load(path+'/ChargePDFs_500r_vhs_fine.npy',allow_pickle=True).tolist()
    e = ImpactEvent([PDFs], np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
    res = e.diffuseCharge(0.0, 0.0, 0.0, 1e8)#, PDFdict['fine'])

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
    
    # v00 = e.getPDF(0, charge_type=hit_type)
    # ax1.contour(gaussian_filter(np.sum(v00.T, axis=0),0.9), norm=mpl.colors.LogNorm(), origin='lower')
    # ax2.contour(gaussian_filter(np.sum(v00.T, axis=1),0.9), norm=mpl.colors.LogNorm(), origin='lower')
    # ax3.contour(gaussian_filter(np.sum(v00.T, axis=2),0.9), norm=mpl.colors.LogNorm(), origin='lower')
    # # ax1.contour(v00.T[:,:,25], norm=mpl.colors.LogNorm(), origin='lower')
    # # ax2.contour(v10.T[:,:,25], norm=mpl.colors.LogNorm(), origin='lower')
    # # ax3.contour(v17.T[:,:,25], norm=mpl.colors.LogNorm(), origin='lower')
    # ax1.set_title('(0,0,0)'); ax2.set_title('(0,0,0.1)'); ax3.set_title('(0,0,0.17)');
    # print e.PDFs['electrons']['dx'],e.PDFs['electrons']['dy'],e.PDFs['electrons']['dz']
if 'pdfs' in run_plots:
    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(halfwidth,halfwidth/2.))
    plot_pdf(ax1,ax2)
    fig.savefig(fig_path+'\pdfs.pdf')

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
        with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        rot_induced = rot_induced[300,0.2]
        plot_rotation_hist(ax_rot, rot_induced, i)
        plot_error(ax_err, rot_induced, i)
    fig.savefig(fig_path+r'\err.pdf')
    
def plot_hist1d(ax, q_induced):
    ax.set_title('Charge Jumps')
    ax.set_xlabel('Jump Size (e)')
    ax.set_ylabel('Counts')
    ax.set_xlim(-0.5,0.5)
    for Q in [1]:
        q = q_induced[300,0.2][:,Q-1]
        r = np.random.normal(0, 0.0132688, q.size)
        e = noiselib.alias(q+r)
        h, bins = np.histogram(e, bins=500, range=(-0.5,0.5))
        x = (bins[1:]+bins[:-1])/2
        center = (bins[:-1] + bins[1:])/2
        center2 = center**2 * np.sign(center)
        bg = np.sum(np.abs(e)>0.05)
        ax.step(center, 6333./1./bg*noiselib.movingmean(h,3), where='mid' )
        # ax.set_ylim([10e-1, 1.5*h.max()])
    ax.set_yscale('log')
if 'hist1d' in run_plots:
    # fig, ax = plt.subplots(1, 1, figsize=(halfwidth.,2), number=700)
    fig = plt.figure(700)
    ax = fig.axes[0]
    for i,hit_type in enumerate(['gammas']):
        with open('dump_sim_impacts_{}_vhs.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        plot_hist1d(ax, q_induced)
    fig.savefig(fig_path+r'\hist1d.pdf')
    
def plot_hist2d(axs, i, q_induced):
    ax1,ax2,ax3 = axs[:,i]
    ax1.set_title(u'$340\ \mathrm{\mu m}$')
    ax2.set_title(u'$640\ \mathrm{\mu m}$')
    ax3.set_title(u'$3195\ \mathrm{\mu m}$');
    # ax2.get_yaxis().set_visible(False); ax3.get_yaxis().set_visible(False);
    ij = [(2,3),(0,1),(0,2)]
    sigma = {'Q1': 0.02086807133980249, 'Q3': 0.008811574164510916, 'Q2': 0.04472867181129812, 'Q4': 0.011038152654827339}
    for j,ax in enumerate([ax1,ax2,ax3]):
        n = q_induced.shape[0]
        r1 = 0.#np.random.normal(0, sigma['Q{}'.format(ij[j][0]+1)], n)
        r2 = 0.#np.random.normal(0, sigma['Q{}'.format(ij[j][1]+1)], n)
        jumps1 = noiselib.alias( q_induced[:,ij[j][0]] + r1, 0.5)
        jumps2 = noiselib.alias( q_induced[:,ij[j][1]] + r2, 0.5)
        ax.scatter(jumps1, jumps2, c='k', marker='.', s=4)
        ax.set_aspect(1)
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,0.5)
        ax.set_xticks([-0.5,-0.25,0,0.25,0.5])
        ax.set_yticks([-0.5,-0.25,0,0.25,0.5])
        ax.set_xticklabels(['','','','',''])
        ax.set_yticklabels(['','','','',''])
        # line_x = [-0.5,-thresh,np.nan,thresh,0.5]
        # line_y = np.full(5,thresh)
        # ax.plot( line_x,  line_y, 'k:' )
        # ax.plot( line_x, -line_y, 'k:' )
        # ax.plot(-line_y,  line_x, 'k:' )
        # ax.plot( line_y,  line_x, 'k:' )
        ax.spines['bottom'].set_linewidth(3)
        ax.spines['left'].set_linewidth(3)
        q1,q2 = [('Q3','Q4'),('Q1','Q2'),('Q1','Q3')][j]
        ax.spines['bottom'].set_color(qcolors[q1])
        ax.spines['left'].set_color(qcolors[q2])
if 'hist2d' in run_plots:
    fig, axs = plt.subplots(3, 2, figsize=(halfwidth,5.5))
    for i,hit_type in enumerate(['gammas','muons']):
        with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
            q_induced, corr, assym, rot_induced = pickle.load(f)
        plot_hist2d(axs, i, q_induced[300,0.2])
    for ax in [axs[2,0]]:
        ax.set_xticklabels([u'\u22120.5','',0,'',0.5])
        ax.set_yticklabels([u'\u22120.5','',0,'',0.5])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xlabel('$\Delta q_\mathrm{1}$ ($e$)')
        ax.set_ylabel('$\Delta q_\mathrm{2}$ ($e$)')
        ax.xaxis.labelpad = 0
        ax.yaxis.labelpad = -12
    fig.savefig(fig_path+r'\hist2d.pdf')

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