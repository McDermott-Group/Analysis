import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys, os
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
import ChargeJumps
reload(ChargeJumps)
from ChargeJumps import *
import itertools

path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
dump_path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data/python_dumps/'
thresh = 0.1
hit_type = 'muons'
L_list = [100,200,300,400,500,600,700,800,900,1000]
fq_list = [ 1.,0.5,0.2,0.1 ]


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
    while pickle.load( open(dump_path+'dump_saveInProgress.dat','rb') ):
        print '.',
        time.sleep(1)
    pickle.dump( True, open(dump_path+'dump_saveInProgress.dat','wb') )
    
def end_file_lock():
    pickle.dump( False, open(dump_path+'dump_saveInProgress.dat','wb') )

def L_fq_iter():
    while True:
        try:
            with open(dump_path+'dump_sim_queue.dat', 'rb') as f:
                L_fq,i = pickle.load(f)
        except IOError:
            L_fq = list(itertools.product(L_list,fq_list))
            i = -1
        i += 1
        with open(dump_path+'dump_sim_queue.dat', 'wb') as f:
            pickle.dump((L_fq,i), f)
        if i >= len(L_fq):
            # os.remove(dump_path+'dump_sim_queue.dat')
            raise StopIteration
        else:
            yield L_fq[i]

def simulate_impacts():
    for L,fq in L_fq_iter():
        
        print('Simulating: {}  L = {}  fQ = {}'.format(hit_type, L, fq))
        
        pdfs_file = load_pdfs(L)
        efiles = load_hits(hit_type)
            
        # run simulation
        app = Controller(sys.argv, fQ=fq, calc_all=True, plot=False,
                         pdfs_file=pdfs_file, event_files=efiles)
        if hasattr(app, 'view'):
            sys.exit(app.exec_())
        while hasattr(app, 'thread') and app.thread.is_alive():
            time.sleep(1)
            
        # save data once sim has finished
        start_file_lock()
        try:
            with open(dump_path+'dump_sim_impacts_{}.dat'.format(hit_type), 'rb') as f:
                data = pickle.load(f)
        except IOError:
            data = { 'q_induced':{}, 'q_direct':{}, 'rot_induced':{},
                'correlation':{}, 'assym1324':{}, 'assym':{}, 'thresh_fraction':{} }
        data['q_induced'][L,fq] =  app.q_induced
        data['q_direct'][L,fq] =  app.q_direct
        data['rot_induced'][L,fq] = app.rot_induced
        with open(dump_path+'dump_sim_impacts_{}.dat'.format(hit_type), 'wb') as f:
            pickle.dump(data, f)
        end_file_lock()
        
        del app
        gc.collect()

def _add_noise(q):
    sigma = [[0.02086807133980249, 0.04472867181129812, 0.008811574164510916, 0.011038152654827339]]
    m_sigma = np.repeat(sigma, q.shape[0], axis=0)
    return q + np.random.normal(0, m_sigma)
    
def add_noise():
    from_file = dump_path+'dump_sim_impacts_{}.dat'.format(hit_type)
    to_file   = dump_path+'dump_sim_impacts_{}_noise.dat'.format(hit_type)
    with open(from_file, 'rb') as f:
        data = pickle.load(f)
    for L in L_list:
        for fq in fq_list:
            data['q_induced'][L,fq] = _add_noise( data['q_induced'][L,fq] )
    with open(to_file, 'wb') as f:
        pickle.dump(data, f)

def gen_qq_plots(fname):
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    for L,fq in data['q_induced'].keys():
        q = data['q_induced'][L,fq]
        fig, axs = plt.subplots(1, 3, figsize=(15,6))
        fig.suptitle('{}:  L = {}  fQ = {}'.format(hit_type, L, fq))
        for i,(q1,q2) in enumerate( ((1,2), (3,4), (1,3)) ):
            ax = axs[i]
            ax.plot( noiselib.alias(q[:,q1-1]),
                     noiselib.alias(q[:,q2-1]), '.' )
            ax.set_xlabel('Q{} aliased charge [e]'.format(q1))
            ax.set_ylabel('Q{} aliased charge [e]'.format(q2))
            ax.set_title('Q{} - Q{}'.format(q1,q2))
            ax.set_aspect(1)
            ax.set_xlim(-0.5,0.5)
            ax.set_ylim(-0.5,0.5)
        plt.draw()
        plt.pause(0.05)
        fig.savefig('{}/qq_figs/{}/L{}fq{}.pdf'.format(path,hit_type,L,fq))
        plt.close(fig)

def calc_corr_assym(fname):
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    CJ = ChargeJumps()
    for L,fq in data['q_induced'].keys():
        q = data['q_induced'][L,fq]
        for i,(q1,q2) in enumerate( ((1,2), (3,4), (1,3)) ):
            data['correlation'][L,fq,(q1,q2)] = \
                    CJ.raw_correlation(q[:,q1-1], q[:,q2-1], thresh, thresh)
            data['assym1324'][L,fq,(q1,q2)] = \
                    CJ.assym1324(q[:,q1-1], q[:,q2-1], thresh, thresh)
        data['assym'][L,fq] = [CJ.assym(q[:,i-1],thresh) for i in (1,2,3,4)]
        data['thresh_fraction'][L,fq] = \
                [CJ.thresh_fraction(q[:,i-1],thresh) for i in (1,2,3,4)]
    with open(fname, 'wb') as f:
        pickle.dump(data, f)


for hit_type in ['gammas']:#,'muons']:
    # simulate_impacts()
    # add_noise()
    # gen_qq_plots(dump_path+'dump_sim_impacts_{}.dat'.format(hit_type))
    # gen_qq_plots(dump_path+'dump_sim_impacts_{}_noise.dat'.format(hit_type))
    # calc_corr_assym( dump_path+'dump_sim_impacts_{}.dat'.format(hit_type) )
    # calc_corr_assym( dump_path+'dump_sim_impacts_{}_noise.dat'.format(hit_type) )
    pass
    
    
""" Print stuff
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
        print( 'Q{} charge asymmetry: {:.3f}'.format( q, 
                            1.*np.sum(e>thresh)/np.sum(np.abs(e)>thresh) ) )
"""

    

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

def plot_dict(d, pair, crange=(0.,1.), measured=None, ax=None, cbo='vertical'):
    
    L_list = np.array([100,200,300,400,500,600,700,800,900,1000])
    fq_list = np.array([1.,0.5,0.2,0.1])
    
    a = np.array([[ d[(L,fq,pair) if pair else (L,fq)] 
                                    for fq in fq_list ] 
                                    for L in L_list ], dtype=float)
    if a.ndim == 3:
        a = np.mean(a, axis=-1)
    
    if measured is not None:
        m = np.nanmax(np.abs(measured - a))
        crange = (measured - m, measured + m)
    p = ax.imshow(a, origin='lower', extent=(-0.5,-0.5+len(fq_list),50,1050), aspect='auto', 
                        vmin=crange[0], vmax=crange[1], cmap='bwr')
    ax.set_ylabel('$\lambda_\mathrm{trap}\ (\mathrm{\mu m})$')
    ax.set_xlabel('$f_q$')
    ax.set_xticks(np.arange(len(fq_list)))
    ax.set_xticklabels(fq_list)
    ar = 10 if cbo == 'horizontal' else 20
    c = ax.figure.colorbar(p, ax=ax, orientation=cbo, aspect=ar)
    c.ax.locator_params(axis='x', nbins=3)
    plt.draw()
    plt.pause(0.05)
        
def del_axes_add_title(axs):
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


with open(dump_path+'dump_sim_impacts_{}_noise.dat'.format('gammas'), 'rb') as f:
    data = pickle.load(f)


""" Plot percentage of events > thresh as a function of L,fQ """
if False:
    
    L0, fq0 = 300, 0.2
    L_list = [100,200,300,400,500,600,700,800,900,1000]
    fq_list = [1,0.5,0.2,0.1]
    
    fig, ax_fq = plt.subplots(1,1)
    ax_fq.set_xlabel('L (um)')
    ax_fq.set_ylabel('Thresh fraction')
    fig, ax_L = plt.subplots(1,1)
    ax_L.set_xlabel('fQ')
    ax_L.set_ylabel('Thresh fraction')
    ax_L.set_xscale('log')

    for ext,c in [('','b'), ('_noise','r')]:

        with open(dump_path+'dump_sim_impacts_{}{}.dat'.format('gammas',ext), 'rb') as f:
            data = pickle.load(f)
            
        d = data['thresh_fraction']
        a = np.array([[ d[(L,fq)] for fq in fq_list ] 
                                  for L in L_list ], dtype=float)
        
        fq_cut = a[:,fq_list.index(fq0),:]
        L_cut = a[L_list.index(L0),:,:]
        ax_fq.plot(L_list, fq_cut, c)
        ax_L.plot(fq_list, L_cut, c)
    
    plt.draw()
    plt.pause(0.05)
