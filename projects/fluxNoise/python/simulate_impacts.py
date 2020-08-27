import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
import matplotlib.pyplot as plt
import noiselib
from numba import jit
import time
import impact_lib
reload(impact_lib)
from impact_lib import *
import pickle

corr = {150:{},300:{},750:{}}
assym = {150:{},300:{},750:{}}
# q_induced = {150:{},300:{},750:{}}

for L in (150, 300, 750):
    for fq in (1., 0.1, 0.01):
        print('Gammas:  L = {}  fQ = {}'.format(L, fq))
        
        app = Controller(sys.argv, fQ=fq, calc_all=False,
                         pdfs_file='sim_data/ChargePDFs_{}.npy'.format(L),
                         event_files=['sim_data/Gamma.txt', 'sim_data/Gamma_10deg.txt'])
                         # event_files=['sim_data/Muons.txt', 'sim_data/Muons2.txt'])
        if hasattr(app, 'view'):
            sys.exit(app.exec_())
        # while app.thread.is_alive():
            # time.sleep(1)
        app.q_induced = q_induced[L][fq]
        print( 'total events: {}'.format(app.q_induced.shape[0]) )
        fig, axs = plt.subplots(1, 3, figsize=(15,6))
        fig.suptitle('Gammas:  L = {}  fQ = {}'.format(L, fq))
        corr[L][fq] = {}
        assym[L][fq] = {}
        for i,(q1,q2) in enumerate( ((1,2), (3,4), (1,3)) ):
            app.plot_qq(q1,q2, axs[i])
            qq1 = app.count(q1,q2,1)
            qq2 = app.count(q1,q2,2)
            qq3 = app.count(q1,q2,3)
            qq4 = app.count(q1,q2,4)
            corr[L][fq][(q1,q2)] = app.get_correlation(q1,q2,0.1)
            try:
                assym[L][fq][(q1,q2)] = 1. * (qq1+qq3) / (qq2+qq4)
            except ZeroDivisionError:
                assym[L][fq][(q1,q2)] = np.nan
            q_induced[L][fq] = app.q_induced
            print( 'Q{} - Q{}'.format(q1,q2) )
            print( '    quadrants 1,2,3,4: {}, {}, {}, {}'.format( qq1, qq2, qq3, qq4 ) )
            try:
                print( '    quadrants 1,2,3,4: {:.2f}, {:.2f}, {:.2f}, {:.2f} %'.format( 
                             *np.array((qq1, qq2, qq3, qq4))/float(qq1+qq2+qq3+qq4) ) )
            except ZeroDivisionError:
                pass
            print( '    correlation: {:.2f}'.format( corr[L][fq][(q1,q2)] ) )
            print( '    13/24 asymmetry: {:.2f}'.format( assym[L][fq][(q1,q2)] ) )
        for q in (1,2,3,4):
            e = noiselib.alias(app.q_induced[:,q-1])
            try:
                print( 'Q{} charge asymmetry: {:.3f}'.format( q, 1.*np.sum(e>0.1)/np.sum(np.abs(e)>0.1) ) )
            except ZeroDivisionError:
                pass
        fig.savefig('sim_data/qq_figs/L{}fq{}.pdf'.format(L,fq))
        plt.close(fig)
        
with open('dump_sim_impacts.dat', 'wb') as f:
    pickle.dump((q_induced,corr,assym), f)

with open('dump_sim_impacts.dat', 'rb') as f:
    q_induced,corr,assym = pickle.load(f)

def print_dict(d):
    print('{:6}{:6}{:10}{:10}{:10}'.format('L','fq','(3,4)','(1,2)','(1,3)'))
    for L in (150, 300, 750):
        for fq in (1., 0.1, 0.01):
            print('{:6}{:6}{:10.2f}{:10.2f}{:10.2f}'.format(
                L, fq, d[L][fq][(3,4)], d[L][fq][(1,2)], d[L][fq][(1,3)] ))
        

fig, axi = plt.subplots(1,1)
axi = plt.figure(1).axes[0]
axi.set_title('Charge Jumps')
axi.set_xlabel('Jump Size [e]')
axi.set_ylabel('')
for l,q_L in q_induced.items():
    for fq,q in q_L.items():
        for Q in [1]:
            e = noiselib.alias(q[:,Q-1])
            h, bins = np.histogram(e, bins=500, range=(-0.5,0.5))
            x = (bins[1:]+bins[:-1])/2
            center = (bins[:-1] + bins[1:])/2
            center2 = center**2 * np.sign(center)
            widths = np.diff(bins**2)
            bg = np.sum(np.abs(e)>0.2)
            axi.step(center2, 2553./bg*noiselib.movingmean(h,30), label='L={} fq={}'.format(l,fq))
axi.set_yscale('log')
axi.set_ylim([10e-1, 1.5*h.max()])
axi.legend()
plt.draw()
plt.pause(0.05)

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
    