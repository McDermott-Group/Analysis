import numpy as np
import matplotlib.pyplot as plt
import noiselib
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
import datasets
reload(datasets)
import datasets as ds


charge_path = '{}\DR1 - 2019-12-17\CorrFar\{}Corr\General\Parameter\{}_correlation.hdf5'
plt.style.use('pub.mplstyle')
fig_path = r'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs'
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



def pubplot_hist_1d():
    CO = make_CO_from_datasets(ds_q1)
    fig, ax = plt.subplots(1, 1, figsize=(halfwidth,2.), num=700)
    jumps, sigma = CO.get_jump_sizes(plot=True, ax=ax, qubits='Q1')
    print np.sum(np.abs(jumps['Q1'])>0.05)
    print sigma

def pubplot_hist_2d():
    CO = make_CO_from_datasets(ds_all4)
    jumps, sigma = CO.get_jump_sizes(plot=False)
    print sigma
    fig,(ax1,ax2,ax3) = plt.subplots(1,3, figsize=(7.2,2.6))
    CO.plot_charge_correlation('Q3','Q4',ax=ax1, thresh=(thresh,thresh))
    CO.plot_charge_correlation('Q1','Q2',ax=ax2, thresh=(thresh,thresh))
    CO.plot_charge_correlation('Q1','Q3',ax=ax3, thresh=(thresh,thresh))
    ax1.set_title(u'$340\ \mathrm{\mu m}$')
    ax2.set_title(u'$640\ \mathrm{\mu m}$')
    ax3.set_title(u'$3195\ \mathrm{\mu m}$');
    # ax2.get_yaxis().set_visible(False); ax3.get_yaxis().set_visible(False);
    for i,ax in enumerate([ax1,ax2,ax3]):
        ax.set_xticks([-0.5,-0.25,0,0.25,0.5])
        ax.set_yticks([-0.5,-0.25,0,0.25,0.5])
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
        ax.spines['bottom'].set_linewidth(3)
        ax.spines['left'].set_linewidth(3)
        q1,q2 = [('Q3','Q4'),('Q1','Q2'),('Q1','Q3')][i]
        ax.spines['bottom'].set_color(qcolors[q1])
        ax.spines['left'].set_color(qcolors[q2])
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


# CO = make_CO_from_datasets(ds_all4)
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_time_steps()

pubplot_hist_1d()
pubplot_hist_2d()

CO = make_CO_from_datasets(ds_all4)
jumps, sigma = CO.get_jump_sizes(plot=False)

""" How many events expected from gamma/muon rates? """
if True:
    rate_muons = 0.5e-3
    rate_gamma = 3e-3
    thresh_fraction_gamma = 0.0690310322989
    thresh_fraction_muons = 0.169750430293
    filter = {q:np.isfinite(jumps[q]) for q in ('Q1','Q2','Q3','Q4')}
    t = CO.plot_time_steps()
    t_q = np.array([np.nansum(t[filter[q]]) for q in ('Q1','Q2','Q3','Q4')])
    print('number of gammas above threshhold:')
    print('    expected = {:.2f}'.format(
        np.nansum(t)*(rate_muons*thresh_fraction_muons 
                    + rate_gamma*thresh_fraction_gamma)))
    # print('    measured = {}'.format(np.sum(np.abs(jumps['Q1'])>thresh)))
    thresh_hits = [np.sum(np.abs(jumps[q])>thresh) for q in ('Q1','Q2','Q3','Q4')]
    rates = (1.*np.array(thresh_hits) - rate_muons*thresh_fraction_muons*t_q)/t_q/1e-3/thresh_fraction_gamma
    print('    measured = {}'.format(thresh_hits))
    print('    measured hit rate = {} mHz'.format(rates))
    print('    mean hit rate = {} mHz'.format(np.mean(rates)))

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
