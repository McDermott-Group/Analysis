import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
from scipy import constants
import matplotlib.pyplot as plt
import noiselib
# from numba import jit
import time
import impact_lib
reload(impact_lib)
from impact_lib import *
import pickle
import gc

path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
calc_all = False
add_noise = False
thresh = 0.1
hit_type = 'gammas'
app = None

def load_pdfs(L):
    # load charge pdf
    if type(L) in (list,tuple):
        pdfs_file = ['{}/ChargePDFs_{}r_vhs_fine.npy'.format(path,l) for l in L]
    else:
        pdfs_file = '{}/ChargePDFs_{}r_vhs_fine.npy'.format(path,L)
    return pdfs_file

def load_hits(hit_type):
    # load generated hits
    if hit_type == 'gammas':
        efiles = [path+'/Gamma.txt', 
                  path+'/Gamma_10deg.txt',
                  path+'/Gamma_10deg_pt2.txt']
    elif hit_type == 'muons':
        efiles = [path+'/Muons.txt', path+'/Muons2.txt']
    return efiles

def start_file_lock():
    while pickle.load( open('dump_saveInProgress.dat','rb') ):
        print '.',
        time.sleep(1)
    pickle.dump( True, open('dump_saveInProgress.dat','wb') )
    
def end_file_lock():
    pickle.dump( False, open('dump_saveInProgress.dat','wb') )

if True:
    for L in [300]:#[100,200,300,400,500,600,700,800,900,1000]:
        for fq in [0.2]:#[ 1.,0.5,0.2,0.1 ]:
        
            print('{}:  L = {}  fQ = {}'.format(hit_type, L, fq))
            
            pdfs_file = load_pdfs(L)
            efiles = load_hits(hit_type)
                
            # run simulation
            app = Controller(sys.argv, fQ=fq, calc_all=calc_all, plot=False,
                             pdfs_file=pdfs_file, event_files=efiles)
            if hasattr(app, 'view'):
                sys.exit(app.exec_())
            if calc_all:
                while hasattr(app, 'thread') and app.thread.is_alive():
                    time.sleep(1)
            
            if not calc_all:
                with open('dump_sim_impacts_{}_noise.dat'.format(hit_type), 'rb') as f:
                    q_induced, corr, assym, rot_induced = pickle.load(f)
                app.q_induced = q_induced[L,fq]
                app.rot_induced = rot_induced[L,fq]
                
            # add noise
            if add_noise:
                sigma = [[0.02086807133980249, 0.04472867181129812, 0.008811574164510916, 0.011038152654827339]]
                m_sigma = np.repeat(sigma, q_induced[L,fq].shape[0], axis=0)
                app.q_induced = app.q_induced + np.random.normal(0, m_sigma)
                
            ### PLOT QQ, CALC ERROR
            start_file_lock()
            try:
                with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
                    q_induced, corr, assym, rot_induced = pickle.load(f)
            except IOError:
                q_induced, corr, assym, rot_induced = {}, {}, {}
            print( 'total events: {}'.format(app.q_induced.shape[0]) )
            fig, axs = plt.subplots(1, 3, figsize=(15,6))
            fig.suptitle('Gammas:  L = {}  fQ = {}'.format(L, fq))
            for i,(q1,q2) in enumerate( ((1,2), (3,4), (1,3)) ):
                app.plot_qq(q1,q2, axs[i])
                qq1 = app.count(q1,q2,1,thresh=thresh)
                qq2 = app.count(q1,q2,2,thresh=thresh)
                qq3 = app.count(q1,q2,3,thresh=thresh)
                qq4 = app.count(q1,q2,4,thresh=thresh)
                corr[L,fq,(q1,q2)] = app.get_correlation(q1,q2,thresh)
                try:
                    a = 1. * (qq1+qq3) / (qq2+qq4)
                except ZeroDivisionError:
                    a = np.nan
                try:
                    da = a * np.sqrt( 1./(qq1+qq3) + 1./(qq2+qq4) )
                except ZeroDivisionError:
                    da = np.nan
                assym[L,fq,(q1,q2)] = a, da
                q_induced[L,fq] = np.array(app.q_induced)
                rot_induced[L,fq] = np.array(app.rot_induced)
                print( 'Q{} - Q{}'.format(q1,q2) )
                print( '    quadrants 1,2,3,4: {}, {}, {}, {}'.format( qq1, qq2, qq3, qq4 ) )
                try:
                    print( '    quadrants 1,2,3,4: {:.2f}, {:.2f}, {:.2f}, {:.2f} %'.format( 
                                 *np.array((qq1, qq2, qq3, qq4))/float(qq1+qq2+qq3+qq4) ) )
                except ZeroDivisionError:
                    pass
                print( u'    correlation: {:.2f} \u00B1 {:.3f}'.format( 
                                *corr[L,fq,(q1,q2)] ) )
                print( u'    13/24 asymmetry: {:.2f} \u00B1 {:.3f}'.format( 
                                *assym[L,fq,(q1,q2)]  ) )
            with open('dump_sim_impacts_{}_noise.dat'.format(hit_type), 'wb') as f:
                pickle.dump((q_induced,corr,assym,rot_induced), f)
            for q in (1,2,3,4):
                e = noiselib.alias(app.q_induced[:,q-1])
                try:
                    print( 'Q{} charge asymmetry: {:.3f}'.format( q, 
                                        1.*np.sum(e>thresh)/np.sum(np.abs(e)>thresh) ) )
                except ZeroDivisionError:
                    pass
            end_file_lock()
            fig.savefig('{}/qq_figs/L{}fq{}.pdf'.format(path,L,fq))
            plt.close(fig)
            del app
            gc.collect()

with open('dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
    q_induced, corr, assym, rot_induced = pickle.load(f)
    

plt.style.use('pub.mplstyle')
fig_path = r'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs'
halfwidth = 3.5
fullwidth = 7.2


""" Plot corr, assym as function of L,fq """
def print_dict(d):
    print(u'{:>6}{:>6}{:>15}{:>15}{:>15}'.format('L','fq','(3,4)','(1,2)','(1,3)'))
    for L in (100, 200, 300, 400, 500, 600, 700, 800):
        for fq in (1., 0.1, 0.01):
            print(u'{:>6}{:>6}{:>15}{:>15}{:>15}'.format(L, fq, 
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(3,4)][0], d[L,fq,(3,4)][1]),
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(1,2)][0], d[L,fq,(1,2)][1]),
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(1,3)][0], d[L,fq,(1,2)][1]) ))

def plot_dict(d, pair, crange=(0.,1.), label='', measured=None, ax=None):
    L_list = np.array([100,200,300,400,500,600,700,800,900,1000])# - 50
    fq_list = np.array([1.,0.5,0.2,0.1])# + 0.5
    data = np.full( (len(L_list),len(fq_list)), 0. )
    for i,L in enumerate(L_list):
        for j,fq in enumerate(fq_list):
            data[i,j] = d[L,fq,pair][0]
    if ax is None:
        fig, ax = plt.subplots(1,1)
    if measured is not None:
        d = np.nanmax(np.abs(measured - data))
        crange = (measured - d, measured + d)
    # p = ax.pcolormesh(L_list, fq_list, data.T, vmin=crange[0], vmax=crange[1])
    p = ax.imshow(data, origin='lower', extent=(-0.5,-0.5+len(fq_list),50,1050), aspect='auto', 
                        vmin=crange[0], vmax=crange[1], cmap='bwr')
    ax.set_ylabel('$\lambda_\mathrm{trap}\ (\mathrm{\mu m})$')
    ax.set_xlabel('$f_q$')
    ax.set_xticks(np.arange(len(fq_list)))
    ax.set_xticklabels(fq_list)
    ax.set_title(label+str(pair))
    c = ax.figure.colorbar(p, ax=ax, orientation='horizontal')
    c.ax.locator_params(axis='x', nbins=3)
    # fig.tight_layout()
    plt.draw()
    plt.pause(0.05)
        
# print_dict(corr)
# print_dict(assym)

# fig, axs = plt.subplots(2,3,constrained_layout=True, figsize=(halfwidth,6.))
# fig.suptitle('Correlation probability ({})'.format(hit_type))
# plot_dict(corr, (3,4), label='Corr ', crange=(0.,0.7), measured=0.54, ax=axs[0,0])
# plot_dict(corr, (1,2), label='Corr ', crange=(0.,0.7), measured=0.46, ax=axs[0,1])
# plot_dict(corr, (1,3), label='Corr ', crange=(0.,0.7), measured=0.0, ax=axs[0,2])
# plot_dict(assym, (3,4), label='13/24 Assym ', crange=(0.,3.5), measured=1.06, ax=axs[1,0])
# plot_dict(assym, (1,2), label='13/24 Assym ', crange=(0.,3.5), measured=1.43, ax=axs[1,1])
# plot_dict(assym, (1,3), label='13/24 Assym ', crange=(0.,3.5), measured=1.27, ax=axs[1,2])
# for ax in np.ravel(axs):
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])
# axs[0,0].set_title('$340\ \mathrm{\mu m}$')
# axs[0,1].set_title('$640\ \mathrm{\mu m}$')
# axs[0,2].set_title('$3195\ \mathrm{\mu m}$')
# axs[1,0].set_title('')
# axs[1,1].set_title('')
# axs[1,2].set_title('')
# plt.draw()
# plt.pause(0.05)
# # fig.savefig(fig_path+'\corr_{}.pdf'.format(hit_type))

# plot correlation probability
if False:
    fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(halfwidth,2.5))
    # fig.suptitle('Correlation probability ({})'.format(hit_type))
    fig.suptitle('Correlation probability')
    plot_dict(corr, (3,4), label='Corr ', crange=(0.,0.7), measured=0.54, ax=axs[0])
    plot_dict(corr, (1,2), label='Corr ', crange=(0.,0.7), measured=0.46, ax=axs[1])
    plot_dict(corr, (1,3), label='Corr ', crange=(0.,0.7), measured=0.0, ax=axs[2])
    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    axs[0].set_title('$340\ \mathrm{\mu m}$')
    axs[1].set_title('$640\ \mathrm{\mu m}$')
    axs[2].set_title('$3195\ \mathrm{\mu m}$')
    plt.draw()
    plt.pause(0.05)
    fig.savefig(fig_path+'\corr_{}.pdf'.format(hit_type))

# plot 13/24 assym
if False:
    fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(halfwidth,2.5))
    # fig.suptitle('13/24 Assymetry ({})'.format(hit_type))
    fig.suptitle('13/24 assymetry')
    plot_dict(assym, (3,4), label='13/24 Assym ', crange=(0.,3.5), measured=1.06, ax=axs[0])
    plot_dict(assym, (1,2), label='13/24 Assym ', crange=(0.,3.5), measured=1.43, ax=axs[1])
    plot_dict(assym, (1,3), label='13/24 Assym ', crange=(0.,3.5), measured=1.27, ax=axs[2])
    for ax in axs:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    axs[0].set_title('$340\ \mathrm{\mu m}$')
    axs[1].set_title('$640\ \mathrm{\mu m}$')
    axs[2].set_title('$3195\ \mathrm{\mu m}$')
    plt.draw()
    plt.pause(0.05)
    fig.savefig(fig_path+r'\assym_{}.pdf'.format(hit_type))

""" Plot percentage of events > thresh as a function of L,fQ """
if True:
    L0, fq0 = 300, 0.2
    L_list = [100,200,300,400,500,600,700,800,900,1000]
    fq_list = [1,0.5,0.2,0.1]
    m_L  = np.array(  [ np.mean(np.abs(q_induced[L, fq0])>thresh,axis=0) for L in L_list])
    m_fq = np.array(  [ np.mean(np.abs(q_induced[L0, fq])>thresh,axis=0) for fq in fq_list])
    m    = np.array( [[ np.mean(np.abs(q_induced[L,  fq])>thresh) 
                        for fq in fq_list] for L in L_list] )
    
    fig, ax = plt.subplots(1,1,constrained_layout=True)
    ax.plot(L_list, m_L)
    ax.set_xlabel('L (um)')
    ax.set_ylabel('% of events > {}$e$'.format(thresh))
    
    fig, ax = plt.subplots(1,1,constrained_layout=True)
    ax.plot(fq_list, m_fq)
    ax.set_xlabel('fQ')
    ax.set_ylabel('% of events > {}$e$'.format(thresh))
    ax.set_xscale('log')
    
    fig, ax = plt.subplots(1,1,constrained_layout=True)
    p = ax.imshow(m, origin='lower', extent=(-0.5,-0.5+len(fq_list),50,1050), aspect='auto')
    fig.colorbar(p, ax=ax)
    ax.set_ylabel('$\lambda_\mathrm{trap}\ (\mathrm{\mu m})$')
    ax.set_xlabel('$f_q$')
    ax.set_xticks(np.arange(len(fq_list)))
    ax.set_xticklabels(fq_list)
    ax.set_title('% events > {}$e$'.format(thresh))
    plt.draw()
    plt.pause(0.05)
    print( '% above thresh for L={},fQ={} = {}'.format(
                L0, fq0, np.mean(np.abs(q_induced[300,0.2])>thresh) ))

""" Overlay 1D histograms """
if False:
    fig, axi = plt.subplots(1,1)
    axi = plt.figure(700).axes[0]
    axi.set_title('Charge Jumps')
    axi.set_xlabel('Jump Size (e)')
    axi.set_ylabel('')
    for (L,fq),q in q_induced.items():
        for Q in [1]:
            e = noiselib.alias(q[:,Q-1])
            h, bins = np.histogram(e, bins=500, range=(-0.5,0.5))
            x = (bins[1:]+bins[:-1])/2
            center = (bins[:-1] + bins[1:])/2
            center2 = center**2 * np.sign(center)
            widths = np.diff(bins**2)
            bg = np.sum(np.abs(e)>0.2)
            axi.step(center2, 2553./bg*noiselib.movingmean(h,30), where='mid',
                        label='L={} fq={}'.format(L,fq))
    axi.set_yscale('log')
    axi.set_ylim([10e-1, 1.5*h.max()])
    noiselib.legend()
    plt.draw()
    plt.pause(0.05)
    
""" Plot +/- assymetry as a function of L,fq """
if True:
    measured = 0.579
    L_list = np.array([100,200,300,400,500,600,700,800,900,1000])# - 50
    fq_list = np.array([1.,0.5,0.2,0.1])# + 0.5
    data = np.full( (len(L_list),len(fq_list),4), 0. )
    for i,L in enumerate(L_list):
        for j,fq in enumerate(fq_list):
            q5 = noiselib.alias(q_induced[L,fq],0.5)
            data[i,j,:] = 1.*np.sum(q5>thresh,axis=0)/np.sum(np.abs(q5)>thresh,axis=0)
    d = np.nanmax(np.abs(measured - data))
    crange = (measured - d, measured + d)
    # fig, axs = plt.subplots(1,4, figsize=(fullwidth,4.))
    # for i in range(4):
        # p = axs[i].imshow(data[:,:,i], origin='lower', 
                          # extent=(-0.5,-0.5+len(fq_list),50,1050), 
                          # vmin=crange[0], vmax=crange[1], aspect='auto', cmap='bwr')
        # axs[i].set_xticks(np.arange(len(fq_list)))
        # axs[i].set_xticklabels(fq_list)
    # fig.colorbar(p, ax=axs)
    # plt.draw()
    # plt.pause(0.05)
    fig, ax = plt.subplots(1,1, figsize=(halfwidth,2.5))
    p = ax.imshow(np.mean(data,axis=2), origin='lower', 
                  extent=(-0.5,-0.5+len(fq_list),50,1050), 
                  aspect='auto', cmap='bwr')
    fig.colorbar(p, ax=ax)
    # ax.set_title('+/- assymetry ({})'.format(hit_type))
    ax.set_title('Charge assymetry')
    ax.set_ylabel('$\lambda_\mathrm{trap}\ \mathrm{(\mu m)}$')
    ax.set_xlabel('$f_q$')
    ax.set_xticks(np.arange(len(fq_list)))
    ax.set_xticklabels(fq_list)
    plt.draw()
    plt.pause(0.05)
    fig.savefig(fig_path+'\pmassym_{}.pdf'.format(hit_type))



# print event.get_induced_charge_on_qubits(np.full((1,3),(0,0,0)), [1])

# if __name__ == '__main__':
    # import numpy as np
    # import matplotlib.pyplot as plt
    # import matplotlib as mpl
    # import simulate_impacts
    # reload(simulate_impacts)
    # from simulate_impacts import ImpactEvent
    # PDFs = np.load('sim_data/ChargePDFs_750.npy',allow_pickle=True).tolist()
    # e = ImpactEvent(PDFs, np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
    # v00 = e.getPDF(0)
    # v10 = e.getPDF(0.1)
    # v17 = e.getPDF(0.17)
    # fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    # ax1.imshow(np.sum(v00.T, axis=1), norm=mpl.colors.LogNorm(), origin='lower')
    # ax2.imshow(np.sum(v10.T, axis=1), norm=mpl.colors.LogNorm(), origin='lower')
    # ax3.imshow(np.sum(v17.T, axis=1), norm=mpl.colors.LogNorm(), origin='lower')
    # ax1.set_title('(0,0,0)'); ax2.set_title('(0,0,0.1)'); ax3.set_title('(0,0,0.17)');
    # plt.draw(); plt.pause(0.05)
