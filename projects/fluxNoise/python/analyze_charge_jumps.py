import numpy as np
import matplotlib.pyplot as plt
import noiselib
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
import ChargeJumps
reload(ChargeJumps)
from ChargeJumps import *
import datasets
reload(datasets)
import datasets as ds
import pickle


charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
plt.style.use('pub.mplstyle')
fig_path = r'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs'
dump_path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data/python_dumps/'
halfwidth = 3.5
fullwidth = 7.2
# qcolors = {'Q1':'C0', 'Q2':'C1', 'Q3':'C2', 'Q4':'C3'}
# qcolors = {'Q1':'firebrick', 'Q2':'lightcoral', 'Q3':'mediumblue', 'Q4':'cornflowerblue'}
qcolors = {'Q1':'blue', 'Q2':'deepskyblue', 'Q3':'red', 'Q4':'magenta'}
thresh = 0.1


ds_all4 = [ds.q1q2q3q4_charge_1, ds.q1q2q3q4_charge_2, ds.q1q2q3q4_charge_3,
           ds.q1q2q3q4_charge_4, ds.q1q2q3q4_charge_5, ds.q1q2q3q4_charge_6,
           ds.q1q2q3q4_charge_7, ds.q1q2q3q4_charge_8, ds.q1q2q3q4_charge_9,
           ds.q1q2q3q4_charge_10]
ds_q1 = [ ds.q1q2_0711_charge_T1, ds.q1q2_0714_charge_T1, ds.q1q2_0715_charge_T1,
          ds.q1q2_0719_charge_T1, ds.q1q2_0720_charge_T1, ds.q1q2_0722_charge_T1,
          ds.q1q2_0722_charge_T1_2, ds.q1q2_0724_charge_T1, ds.q1q2_0724_charge_T1_2,
          ds.q1q2_0725_charge_T1, ds.q1q2q4_0726_charge_T1, ds.q1q2q4_0727_charge_T1, 
          ds.q1q2q4_0728_charge_T1, ds.q1q2q4_0729_charge_T1, ds.q1q2q4_0730_charge_T1, 
          ds.q1q2q4_0731_charge_T1, ds.q1q2q4_0801_charge_T1, ds.q1q2q4_0802_charge_T1, 
          # ds.q1q4_0803_charge_T1, ds.q1q4_0805_charge_T1, 
          ds.q1q2q4_0810_charge_T1, ds.q1q2q4_0811_charge_T1, ds.q1q2q4_0812_charge_T1] # Q1 only
ds_all4_2 = [ds.q1q2q3q4_charge_11, ds.q1q2q3q4_charge_12, ds.q1q2q3q4_charge_13, 
             ds.q1q2q3q4_charge_14]

def make_CO_from_datasets(ds_list):
    CO = ChargeOffset()
    for d in ds_list:
        if 'path_base' in d:
            base = d['path_base']
        else:
            base = 'fluxNoise2'
        CO.add_dataset(charge_path.format(base, d['Q'], d['path_charge']))
        if 'exclude_data' in d and type(d['exclude_data'])==dict:
            for q,ex in d['exclude_data'].items():
                for start,end in ex:
                    CO.limit_dataset(q,start=start,end=end)
        if 'exclude_data' in d and type(d['exclude_data'])==list:
            for start,end in d['exclude_data']:
                CO.limit_dataset('Q1',start=start,end=end)
            try:
                CO.limit_dataset('Q2')
                CO.limit_dataset('Q3')
                CO.limit_dataset('Q4')
            except KeyError:
                pass
    return CO

def save_datasets():
    CO = make_CO_from_datasets(ds_q1)
    jumps, sigma = CO.get_jump_sizes(plot=True, ax=ax, qubits='Q1')
    q = jumps['Q1']
    with open(dump_path+'dump_measured_Q1.dat', 'wb') as f:
        pickle.dump(q, f)
    CO = make_CO_from_datasets(ds_all4)
    jumps, sigma = CO.get_jump_sizes(plot=False)
    q = np.vstack([jumps['Q1'],jumps['Q2'],jumps['Q3'],jumps['Q4']]).T
    with open(dump_path+'dump_measured_Q1234.dat', 'wb') as f:
        pickle.dump(q, f)


def pubplot_err_1d():
    CO = make_CO_from_datasets(ds_q1)
    fig, ax = plt.subplots(1, 1, figsize=(fullwidth,3.), num=700)
    jumps, sigma = CO.get_jump_sizes(qubits='Q1')
    dw01 = 2. * np.pi * 12e3
    Ec = 2. * np.pi * 250e6
    tau = 1e-6
    dq = jumps['Q1'][ np.isfinite(jumps['Q1']) ]
    err = dw01**2 / 12. * np.sin(np.pi*dq/2.)**2 * tau**2
    
    h, bins = np.histogram(err, bins=np.logspace(np.log10(1e-14),np.log10(1e-3),101))
    ax.bar( bins[:-1], 1.*h/np.sum(h), width=np.diff(bins), align='edge' )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$\epsilon_\phi$')
    ax.set_ylabel('Normalized Counts')
    plt.draw()
    plt.pause(0.05)

def pubplot_err_2d():
    fig,(ax1,ax2,ax3) = plt.subplots(1,3, figsize=(fullwidth,2.6))
    CO = make_CO_from_datasets(ds_all4)
    jumps, sigma = CO.get_jump_sizes(plot=False)
    dw01 = 2. * np.pi * 12e3
    Ec = 2. * np.pi * 250e6
    tau = 1e-6
    for i,ax in enumerate((ax1,ax2,ax3)):
        j,k = [(3,4),(1,2),(1,3)][i]
        dq1, dq2 = noiselib.overlap( jumps['Q%i'%j], jumps['Q%i'%k] )
        err1 = dw01**2 / 12. * np.sin(np.pi*dq1/2.)**2 * tau**2
        err2 = dw01**2 / 12. * np.sin(np.pi*dq2/2.)**2 * tau**2
        
        ax.scatter(err1[np.abs(dq1)<0.1], err2[np.abs(dq1)<0.1], marker='.', s=4)
        ax.scatter(err1[np.abs(dq1)>0.1], err2[np.abs(dq1)>0.1], marker='.', s=4)
        ax.scatter(err1[np.abs(dq2)>0.1], err2[np.abs(dq2)>0.1], marker='.', s=4)
        # ax.scatter(err1, err2, c='k', marker='.', s=4)
        noiselib.format_hist2d(i, ax, log=True, cprofile=qcolors)
    noiselib.edges_tick_labels([ax1,ax2,ax3], xlabel=u'$\epsilon_1$', ylabel=u'$\epsilon_1$')
    
    plt.draw()
    plt.pause(0.05)

def pubplot_hist_1d():
    CO = make_CO_from_datasets(ds_q1)
    fig, ax = plt.subplots(1, 1, figsize=(fullwidth,3.), num=700)
    jumps, sigma = CO.get_jump_sizes(plot=True, ax=ax, qubits='Q1')
    ax.lines[1].set_visible(False)
    print np.sum(np.abs(jumps['Q1'])>0.05)
    print sigma

def pubplot_hist_2d():
    CO = make_CO_from_datasets(ds_all4)
    jumps, sigma = CO.get_jump_sizes(plot=False)
    print sigma
    fig,(ax1,ax2,ax3) = plt.subplots(1,3, figsize=(fullwidth,2.6))
    CJ = ChargeJumps()
    print CO.plot_charge_correlation('Q3','Q4',ax=ax1, thresh=(thresh,thresh), plot=False)
    print CJ.plot_charge_correlation(CO, 'Q3','Q4',ax=ax1, thresh=(thresh,thresh), plot=False)
    _,_,_,_,g3,dg3,g4,dg4 = CO.plot_charge_correlation('Q3','Q4',ax=ax1, thresh=(thresh,thresh))
    _,_,_,_,g1,dg1,g2,dg2 = CO.plot_charge_correlation('Q1','Q2',ax=ax2, thresh=(thresh,thresh))
    CO.plot_charge_correlation('Q1','Q3',ax=ax3, thresh=(thresh,thresh))
    # ax2.get_yaxis().set_visible(False); ax3.get_yaxis().set_visible(False);
    for i,ax in enumerate([ax1,ax2,ax3]):
        noiselib.format_hist2d(i, ax, cprofile=qcolors)
        ax.set_xticklabels([u'\u22120.5','',0,'',0.5])
        ax.set_yticklabels(['','','','',''])
        ax.set_xlabel('')
        ax.set_ylabel('')
        line_x = [-0.5,-thresh,np.nan,thresh,0.5]
        line_y = np.full(5,thresh)
        ax.plot( line_x,  line_y, 'k:' )
        ax.plot( line_x, -line_y, 'k:' )
        ax.plot(-line_y,  line_x, 'k:' )
        ax.plot( line_y,  line_x, 'k:' )
    ax1.set_yticklabels([u'\u22120.5','',0,'',0.5])
    ax1.set_xlabel('$\Delta q_\mathrm{1}$ ($e$)')
    ax1.set_ylabel('$\Delta q_\mathrm{2}$ ($e$)')
    ax1.xaxis.labelpad = 0
    ax1.yaxis.labelpad = -12
    fig.savefig(fig_path+'\qq.pdf')
    print 'jump assymetry:'
    for q in ('Q1','Q2','Q3','Q4'):
        print '    {}: {:.3f}'.format( q, 1.*np.sum(jumps[q]>thresh)/np.sum(np.abs(jumps[q])>thresh) )
    plt.draw()
    plt.pause(0.05)
    return np.array([g1,g2,g3,g4]), np.array([dg1,dg2,dg3,dg4])


# CO = make_CO_from_datasets(ds_all4)
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_time_steps()

# save_datasets()
# pubplot_err_1d()
# pubplot_err_2d()
# pubplot_hist_1d()
# rates, d_rates = pubplot_hist_2d()

CO = make_CO_from_datasets(ds_all4)
jumps, sigma = CO.get_jump_sizes(plot=False)
dt = CO.get_time_steps()
CJ = ChargeJumps()

p_obs, p, G, pAB_norm = {}, {}, {}, {}
for q1,q2 in [('Q3','Q4'),('Q1','Q2')]:
    (p_obs[q1], p_obs[q2], p[q1], p[q2], G[q1], G[q2], 
     p_obs[(q1,q2)], p[(q1,q2)], G[(q1,q2)], pAB_norm[(q1,q2)], a1324) \
     = CJ.calc_params(jumps[q1], jumps[q2], thresh, thresh, dt)
for q1,q2 in [('Q1','Q3')]:
    ( _,_,_,_,_,_, 
     p_obs[(q1,q2)], p[(q1,q2)], G[(q1,q2)], pAB_norm[(q1,q2)], a1324) \
     = CJ.calc_params(jumps[q1], jumps[q2], thresh, thresh, dt)
     
for k in ('Q1','Q2','Q3','Q4',('Q3','Q4'),('Q1','Q2'),('Q1','Q3')):
    print( '{:>15} {:>15} {:>15} {:>15}'.format(k, p_obs[k], p[k], G[k]) )
for k in (('Q3','Q4'),('Q1','Q2'),('Q1','Q3')):
    print( '{:>15} {:>15}'.format(k, pAB_norm[k]) )

def get_thresh_fractions():
    L0, fq0 = 300, 0.2
    data_files = {}
    data_files['gammas'] = [dump_path+'dump_sim_impacts_{}_noise_{}.dat'.format('gammas',i)
                                        for i in (0,1,2,3)]
    data_files['muons'] = [dump_path+'dump_sim_impacts_{}_noise_0.dat'.format('muons')]
    thresh_frac = {}
    for h,hit_type in enumerate(['gammas', 'muons']):
        data = []
        for fname in data_files[hit_type]:
            with open(fname, 'rb') as f:
                data += [pickle.load(f)]
        thresh_frac[hit_type] = np.mean( [ d['thresh_fraction'][L0,fq0]
                                                for d in data ] )
    return thresh_frac

""" How many events expected from gamma/muon rates? """
if True:
    rate_muons = 0.5e-3
    rate_gamma = 3e-3
    thresh_frac = get_thresh_fractions()
    rates = np.array([G[q] for q in ('Q1','Q2','Q3','Q4')])
    d_rates = np.array([G[q].e for q in ('Q1','Q2','Q3','Q4')])
    avg_rate = np.sum(rates/d_rates**2) / np.sum(1/d_rates**2)
    avg_rate_gammas = avg_rate - rate_muons*thresh_frac['muons']
    print( rates/1e-3 )
    print( 'Simulated thresh fractions: {}'.format(thresh_frac) )
    print( 'Avg rate above thresh (all): {} mHz'.format(avg_rate/1e-3) )
    print( 'Avg rate above thresh (gammas): {} mHz'.format(avg_rate_gammas/1e-3) )
    print( 'Avg total rate (gammas): {} mHz'.format(
                avg_rate_gammas/thresh_frac['gammas']/1e-3) )

# thresh_list = np.arange(0.05, 0.35, 0.01)
# n = thresh_list.size
# pC = { 12: np.full((n,n),np.nan),
       # 34: np.full((n,n),np.nan),
       # 13: np.full((n,n),np.nan) }
# d_pC = { 12: np.full((n,n),np.nan),
       # 34: np.full((n,n),np.nan),
       # 13: np.full((n,n),np.nan) }
# a1324 = { 12: np.full((n,n),np.nan),
       # 34: np.full((n,n),np.nan),
       # 13: np.full((n,n),np.nan) }
# d_a1324 = { 12: np.full((n,n),np.nan),
       # 34: np.full((n,n),np.nan),
       # 13: np.full((n,n),np.nan) }
# for i, thresh1 in enumerate(thresh_list):
    # for j, thresh2 in enumerate(thresh_list):
        # print i,j
        # pC[12][i,j], d_pC[12][i,j], a1324[12][i,j], d_a1324[12][i,j] = \
            # CO.plot_charge_correlation('Q1','Q2',ax=ax1, thresh=(thresh1,thresh2), plot=False)
        # pC[34][i,j], d_pC[34][i,j], a1324[34][i,j], d_a1324[34][i,j] = \
            # CO.plot_charge_correlation('Q3','Q4',ax=ax1, thresh=(thresh1,thresh2), plot=False)
        # pC[13][i,j], d_pC[13][i,j], a1324[13][i,j], d_a1324[13][i,j] = \
            # CO.plot_charge_correlation('Q1','Q3',ax=ax1, thresh=(thresh1,thresh2), plot=False)


# for qq in [34,12,13]:
    # fig,axs = plt.subplots(1, 4, figsize=(15,4.5))
    # fig.suptitle(str(qq))
    # for i,(name,var) in enumerate([('pC',pC), ('error pC',d_pC), 
                                   # ('13/24',a1324), ('error 13/24',d_a1324)]):
        # m = axs[i].pcolormesh(thresh_list, thresh_list, var[qq])
        # axs[i].set_aspect( 1 )
        # axs[i].set_title( name )
        # plt.colorbar(m, ax=axs[i], location='bottom')
        # plt.draw(); plt.pause(0.05)

# thresh_list = np.arange(0.05, 0.35, 0.01)
# plt.rcParams['errorbar.capsize'] = 2
# fig1, ax1 = plt.subplots(1, 1)
# fig2, ax2 = plt.subplots(1, 1)
# for qq in [34,12,13]:
    # ax1.errorbar( thresh_list, np.diag(pC[qq]), np.diag(d_pC[qq]), fmt='.', label=str(qq))
    # ax2.errorbar( thresh_list, np.diag(a1324[qq]), np.diag(d_a1324[qq]), fmt='.', label=str(qq))
# ax1.legend()
# ax2.legend()
# ax1.set_title('pC')
# ax2.set_title('13/24')
# plt.draw()
# plt.pause(0.05)

# jumps, sigma = CO.get_jump_sizes()
# fig, axs = plt.subplots(1, 3)
# for i,(q1,q2) in enumerate([('Q3','Q4'),('Q1','Q2'),('Q1','Q3')]):
    # corrJumps = [ np.sum(CO._above(jumps[q1], t) & CO._above(jumps[q2], t)) for t in thresh_list]
    # jumps1 = [ np.sum(CO._above(jumps[q1], t) & CO._above(jumps[q2], 0)) for t in thresh_list]
    # jumps2 = [ np.sum(CO._above(jumps[q2], t) & CO._above(jumps[q1], 0)) for t in thresh_list]
    # axs[i].plot( thresh_list, jumps1, label='N Jumps {}'.format(q1) )
    # axs[i].plot( thresh_list, jumps2, label='N Jumps {}'.format(q2) )
    # axs[i].plot( thresh_list, corrJumps, label='N Corr Jumps' )
    # axs[i].legend()
    # axs[i].set_title('{} - {}'.format(q1,q2))
    # axs[i].set_xlabel('Threshhold')
# plt.draw()
# plt.pause(0.05)

"""######### Interleaved Data ##########"""

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q3Q4Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csn0100hzy_correlation.hdf5')
# CO.add_dataset(base_path + 'csn1803vvg_correlation.hdf5')
# CO.add_dataset(base_path + 'cso0030jjg_correlation.hdf5')
# CO.add_dataset(base_path + 'csq0128zge_correlation.hdf5')
# # CO.add_dataset(base_path + 'csy2025oeg_correlation.hdf5')
# # CO.add_dataset(base_path + 'csz0145znf_correlation.hdf5')
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q3','Q4')
# CO.plot_time_steps()
# plt.draw()
# plt.pause(0.05)

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csu0052mrt_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0034shs_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0824hxv_correlation.hdf5')
# # CO.add_dataset(base_path + 'cvh0508hpk_correlation.hdf5') # dif timing, etc?  good data though
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q2')
# CO.plot_time_steps()
# plt.draw()
# plt.pause(0.05)

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q3Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csw0353opr_correlation.hdf5')
# CO.add_dataset(base_path + 'csx0545cpf_correlation.hdf5')
# CO.add_dataset(base_path + 'csx1710hcz_correlation.hdf5')
# CO.add_dataset(base_path + 'csy0726yvn_correlation.hdf5')
# # base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q4\General\Parameter\\'
# # CO.add_dataset(base_path + 'crc2345hmc_parameters.hdf5') # might say Q3,Q4 instead?
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q3')
# CO.plot_time_steps()
# plt.draw()
# plt.pause(0.05)
