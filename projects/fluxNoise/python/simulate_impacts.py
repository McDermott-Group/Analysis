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
import gc

path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
calc_all = False
thresh = 0.12


# # for L in [100, 200, 300, 400, 500, 600, 700, 800]:
# for L in [800]:
    # for fq in [ 0.01 ]:
    
        # ### SIMULATE
        # print('Gammas:  L = {}  fQ = {}'.format(L, fq))
        
        # if type(L) in (list,tuple):
            # pdfs_file = ['{}/ChargePDFs_{}.npy'.format(path,l) for l in L]
        # else:
            # pdfs_file = '{}/ChargePDFs_{}.npy'.format(path,L)
        # app = Controller(sys.argv, fQ=fq, calc_all=calc_all, plot=False,
                         # pdfs_file=pdfs_file,
                         # event_files=[path+'/Gamma.txt', 
                                      # path+'/Gamma_10deg.txt',
                                      # path+'/Gamma_10deg_pt2.txt'])
                         # # event_files=[path+'/Muons.txt', path+'sim_data/Muons2.txt'])
        # if hasattr(app, 'view'):
            # sys.exit(app.exec_())
        # if not calc_all:
            # with open('dump_sim_impacts.dat', 'rb') as f:
                # q_induced, corr, assym = pickle.load(f)
            # app.q_induced = q_induced[L,fq]
        # else:
            # while hasattr(app, 'thread') and app.thread.is_alive():
                # time.sleep(1)
            
        # ### PLOT QQ, CALC ERROR
        # while pickle.load( open('dump_saveInProgress.dat','rb') ):
            # print '.',
            # time.sleep(1)
        # pickle.dump( True, open('dump_saveInProgress.dat','wb') )
        # with open('dump_sim_impacts.dat', 'rb') as f:
            # q_induced, corr, assym = pickle.load(f)
        # print( 'total events: {}'.format(app.q_induced.shape[0]) )
        # fig, axs = plt.subplots(1, 3, figsize=(15,6))
        # fig.suptitle('Gammas:  L = {}  fQ = {}'.format(L, fq))
        # for i,(q1,q2) in enumerate( ((1,2), (3,4), (1,3)) ):
            # app.plot_qq(q1,q2, axs[i])
            # qq1 = app.count(q1,q2,1,thresh=thresh)
            # qq2 = app.count(q1,q2,2,thresh=thresh)
            # qq3 = app.count(q1,q2,3,thresh=thresh)
            # qq4 = app.count(q1,q2,4,thresh=thresh)
            # corr[L,fq,(q1,q2)] = app.get_correlation(q1,q2,thresh)
            # try:
                # a = 1. * (qq1+qq3) / (qq2+qq4)
            # except ZeroDivisionError:
                # a = np.nan
            # try:
                # da = a * np.sqrt( 1./(qq1+qq3) + 1./(qq2+qq4) )
            # except ZeroDivisionError:
                # da = np.nan
            # assym[L,fq,(q1,q2)] = a, da
            # q_induced[L,fq] = np.array(app.q_induced)
            # print( 'Q{} - Q{}'.format(q1,q2) )
            # print( '    quadrants 1,2,3,4: {}, {}, {}, {}'.format( qq1, qq2, qq3, qq4 ) )
            # try:
                # print( '    quadrants 1,2,3,4: {:.2f}, {:.2f}, {:.2f}, {:.2f} %'.format( 
                             # *np.array((qq1, qq2, qq3, qq4))/float(qq1+qq2+qq3+qq4) ) )
            # except ZeroDivisionError:
                # pass
            # print( u'    correlation: {:.2f} \u00B1 {:.3f}'.format( 
                            # *corr[L,fq,(q1,q2)] ) )
            # print( u'    13/24 asymmetry: {:.2f} \u00B1 {:.3f}'.format( 
                            # *assym[L,fq,(q1,q2)]  ) )
        # with open('dump_sim_impacts.dat', 'wb') as f:
            # pickle.dump((q_induced,corr,assym), f)
        # for q in (1,2,3,4):
            # e = noiselib.alias(app.q_induced[:,q-1])
            # try:
                # print( 'Q{} charge asymmetry: {:.3f}'.format( q, 
                                    # 1.*np.sum(e>thresh)/np.sum(np.abs(e)>thresh) ) )
            # except ZeroDivisionError:
                # pass
        # pickle.dump( False, open('dump_saveInProgress.dat','wb') )
        # fig.savefig('{}/qq_figs/L{}fq{}.pdf'.format(path,L,fq))
        # plt.close(fig)
        # del app
        # gc.collect()

with open('dump_sim_impacts.dat', 'rb') as f:
    q_induced, corr, assym = pickle.load(f)
    
    
""" Plot +/- assymetry as a function of L,fq """
data = np.full( (8,3,4), 0. )
for i,L in enumerate((100, 200, 300, 400, 500, 600, 700, 800)):
    for j,fq in enumerate((1., 0.1, 0.01)):
        q5 = noiselib.alias(q_induced[L,fq],0.5)
        data[i,j,:] = 1.*np.sum(q5>thresh,axis=0)/np.sum(np.abs(q5)>thresh,axis=0)
fig, axs = plt.subplots(1,4,constrained_layout=True)
for i in range(4):
    p = axs[i].imshow(data[:,:,i], origin='lower', extent=(0.5,-2.5,50,850), 
                      vmin=0.3, vmax = 0.7, aspect='auto', cmap='bwr')
fig.colorbar(p, ax=axs)
plt.draw()
plt.pause(0.05)
fig, ax = plt.subplots(1,1,constrained_layout=True)
p = ax.imshow(np.mean(data,axis=2), origin='lower', extent=(0.5,-2.5,50,850), 
                  aspect='auto', cmap='bwr')
fig.colorbar(p, ax=ax)
ax.set_ylabel('L')
ax.set_xlabel('fq $10^x$')
ax.set_title('Avg +/- assympetry')
plt.draw()
plt.pause(0.05)


""" Plot corr, assym as function of L,fq """
def print_dict(d):
    print(u'{:>6}{:>6}{:>15}{:>15}{:>15}'.format('L','fq','(3,4)','(1,2)','(1,3)'))
    for L in (100, 200, 300, 400, 500, 600, 700, 800):
        for fq in (1., 0.1, 0.01):
            print(u'{:>6}{:>6}{:>15}{:>15}{:>15}'.format(L, fq, 
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(3,4)][0], d[L,fq,(3,4)][1]),
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(1,2)][0], d[L,fq,(1,2)][1]),
                u'{:.2f} \u00B1 {:.3f}'.format(d[L,fq,(1,3)][0], d[L,fq,(1,2)][1]) ))

def plot_dict(d, pair, range=(0.,1.), label='', measured=None):
    L_list = np.array([100,200,300,400,500,600,700,800])# - 50
    fq_list = np.array([0,-1,-2])# + 0.5
    data = np.full( (len(L_list),len(fq_list)), 0. )
    for i,L in enumerate([100,200,300,400,500,600,700,800]):
        for j,fq in enumerate([1.,0.1,0.01]):
            data[i,j] = d[L,fq,pair][0]
    fig, ax = plt.subplots(1,1,constrained_layout=True)
    if measured is not None:
        d = np.nanmax(np.abs(measured - data))
        range = (measured - d, measured + d)
    # p = ax.pcolormesh(L_list, fq_list, data.T, vmin=range[0], vmax=range[1])
    p = ax.imshow(data, origin='lower', extent=(0.5,-2.5,50,850), aspect='auto', 
                        vmin=range[0], vmax=range[1], cmap='bwr')
    ax.set_ylabel('L')
    ax.set_xlabel('fq $10^x$')
    ax.set_title(label+str(pair))
    fig.colorbar(p, ax=ax)
    # fig.tight_layout()
    plt.draw()
    plt.pause(0.05)
    return L_list, fq_list, data
        
print_dict(corr)
print_dict(assym)

plot_dict(corr, (3,4), label='Corr ', range=(0.,0.7), measured=0.52)
plot_dict(corr, (1,2), label='Corr ', range=(0.,0.7), measured=0.43)
plot_dict(corr, (1,3), label='Corr ', range=(0.,0.7), measured=0.0)
plot_dict(assym, (3,4), label='13/24 Assym ', range=(0.,3.5), measured=1.03)
plot_dict(assym, (1,2), label='13/24 Assym ', range=(0.,3.5), measured=1.38)
plot_dict(assym, (1,3), label='13/24 Assym ', range=(0.,3.5), measured=1.5)

""" Overlay 1D histograms """
# fig, axi = plt.subplots(1,1)
# axi = plt.figure(1).axes[0]
# axi.set_title('Charge Jumps')
# axi.set_xlabel('Jump Size [e]')
# axi.set_ylabel('')
# for (L,fq),q in q_induced.items():
    # for Q in [1]:
        # e = noiselib.alias(q[:,Q-1])
        # h, bins = np.histogram(e, bins=500, range=(-0.5,0.5))
        # x = (bins[1:]+bins[:-1])/2
        # center = (bins[:-1] + bins[1:])/2
        # center2 = center**2 * np.sign(center)
        # widths = np.diff(bins**2)
        # bg = np.sum(np.abs(e)>0.2)
        # axi.step(center2, 2553./bg*noiselib.movingmean(h,30),
                    # label='L={} fq={}'.format(L,fq))
# axi.set_yscale('log')
# axi.set_ylim([10e-1, 1.5*h.max()])
# noiselib.legend()
# plt.draw()
# plt.pause(0.05)

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
