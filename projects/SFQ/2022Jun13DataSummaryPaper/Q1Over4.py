# from SFQlib import RB, RB_AllGates, Purity, T1_QP_2D_Linear, T1_QP_2D, add_2Ddata_from_matlab
# from SFQlib import RB_AllGates_Paper, Purity_Paper, Purity_Paper_ErrorBudget, T1_QP_Paper
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

if 0:  # IRB the optimized result 1.19% error/clifford gate
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022May30RB'  #
    exp_name = 'RB_All'
    file_Number = np.arange(0, 6, 1)
    RB_file = [
        file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i)
        for i in file_Number]
    RB_data = RB_AllGates_Paper()
    RB_data.add_data_from_matlab(RB_file)
    RB_data.data_analysis()
    RB_data.plot()

if 0:  # purity data, with 0.964 % incoherence error/clifford
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022Jun05_PurityVsRB'
    exp_name = 'Purity2D'
    file_Number = np.arange(0, 10, 1)
    file = [
        file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i)
        for i in file_Number]
    # Purity_data = Purity_Paper()
    Purity_data = Purity_Paper_ErrorBudget()
    Purity_data.add_data_from_matlab(file)
    Purity_data.data_analysis()
    Purity_data.plot()

if 0:  # compare the fidelity data
    data = np.array([
        [10, 99.867, 0.084], [20, 99.556, 0.099], [30, 99.552, 0.121],
        [40, 99.364, 0.178], [50, 99.427, 0.156], [60, 99.292, 0.140],
        [70, 99.118, 0.116], [80, 98.823, 0.185]])

    T1_base = 26.0  # us
    Twhite_base = 2 * T1_base  # us
    Texp_base = 3 / (1 / T1_base + 1 / Twhite_base)
    print('Texp_base=', Texp_base)
    k_base = 1 / (Texp_base * 1000) * 100  # from Zijun Chen's thesis

    r_incoherent = 0.964  # 0.964% per clifford gate from purity
    r_RB = 1.19  # 0.964% per clifford gate from RB/IRB

    time = data[:, 0]
    r = 100 - data[:, 1]
    se = data[:, 2]
    [[k, b], cov] = np.polyfit(time, r, 1, cov=True)
    std = np.sqrt(np.diag(cov))

    plt.errorbar(time, r, yerr=se, fmt="o", color='k', capsize=3.0, ms=5.0)
    plt.plot(time, time * k + b, 'k--',
             label='Idle gate error= {:.3f}  $\pm$ {:.3f} % (per 10 ns)'
             .format(k * 10, std[0] * 10))

    # plt.plot(time, time * r_RB/91, 'y--', label='IRB error= {:.3f} % (per 10 ns)'
    #          .format(r_RB/91 * 10))

    plt.plot(time, time * r_incoherent / 91, 'r--',
             label='Purity incoherent error= {:.3f} % (per 10 ns)'
             .format(r_incoherent / 91 * 10))

    plt.plot(time, time * k_base, 'b--',
             label='Base incoherent error= {:.3f} % (per 10 ns)'
             .format(k_base * 10))
    # plt.errorbar(79, 100-99.151, yerr=0.169, fmt="s", color='r', capsize=3.0, ms=8, label='X')
    # plt.errorbar(79, 100-99.125, yerr=0.179, fmt="v", color='b', capsize=3.0, ms=8, label='Y')
    # plt.errorbar(39, 100-99.503, yerr=0.160, fmt="*", color='r', capsize=3.0, ms=12, label='X/2')
    # plt.errorbar(39, 100-99.383, yerr=0.128, fmt="D", color='b', capsize=3.0, ms=8, label='-X/2')
    # plt.errorbar(39, 100-99.377, yerr=0.180, fmt="p", color='r', capsize=3.0, ms=10, label='Y/2')
    # plt.errorbar(39, 100-99.338, yerr=0.169, fmt="h", color='b', capsize=3.0, ms=10, label='-Y/2')

    plt.plot(79, 100 - 99.151, marker="s", color='r', ms=7, label='X')
    plt.plot(39, 100 - 99.503, marker="*", color='r', ms=10, label='X/2')
    plt.plot(39, 100 - 99.383, marker="D", color='r', ms=6, label='-X/2')
    plt.plot(79, 100 - 99.125, marker="v", color='b', ms=7, label='Y')
    plt.plot(39, 100 - 99.377, marker="p", color='b', ms=8, label='Y/2')
    plt.plot(39, 100 - 99.338, marker="h", color='b', ms=8, label='-Y/2')

    plt.title('Identity Gate Error vs Gate Time')

    plt.legend()
    plt.xlabel('Identity Gate Time (ns)')
    plt.ylabel('Identity Gate Error (%)')
    plt.show()

if 0:  # Q1 T1 with SFQ1 off-resonant poison
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022Jun01_QPStudy'
    experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep'
    file_Number = np.arange(1, 13, 1)
    T1_2D_file = [file_path.format(date, experiment_name_T1,
                                   experiment_name_T1) + '_{:03d}.mat'.format(
        i) for i in file_Number]
    T1_2D_data = T1_QP_2D_Linear()
    T1_2D_data.add_data_from_matlab(T1_2D_file)
    T1_2D_data.plot()

if 0:  # Q1 T1 data with SFQ2 Poisoning
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '06-02-22'
    experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep_Q1'
    file_Number = [0, 1, 2, 3, 4]
    T1_2D_file = [file_path.format(date, experiment_name_T1,
                                   experiment_name_T1) + '_{:03d}.mat'.format(
        i) for i in file_Number]
    T1_2D_data = T1_QP_2D_Linear()
    T1_2D_data.add_data_from_matlab(T1_2D_file)
    T1_2D_data.plot()

if 0:  # Q1 T1 data with SFQ2 Poisoning finer steps
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = '2022Jun01_QPStudy'
    experiment_name_T1 = 'T1_Q1_SFQ2_Poison'
    file_Number = [0, 1, 2, 3, 4, 9]
    T1_2D_file = [file_path.format(date, experiment_name_T1,
                                   experiment_name_T1) + '_{:03d}.mat'.format(
        i) for i in file_Number]
    T1_2D_data = T1_QP_2D_Linear()
    T1_2D_data.add_data_from_matlab(T1_2D_file)
    T1_2D_data.plot()

if 0:  # QP data
    #        p0 = [0.2, 6e3, 20e3, 0.95, 0.1]  # initial guess
    #     bounds = [(0, 5e3, 15e3, 0.5, 0), (2, 9e3, 30e3, 1.05, 0.3)]  # bounds
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = 'T1PoisonSweep'
    experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep_Over4'
    n_QP_2D = np.array([])
    for i in [0, 1, 2, 3, 4, 5, 6, 7]:
        file_Number = [i]
        T1_2D_file = [file_path.format(date, experiment_name_T1,
                                       experiment_name_T1) + '_{:03d}.mat'.format(
            i) for i in file_Number]
        T1_2D_data = T1_QP_2D()
        T1_2D_data.add_data_from_matlab(T1_2D_file)
        n_QP_1D = T1_2D_data.params_2D[:, 0]
        if len(n_QP_2D) == 0:
            n_QP_2D = n_QP_1D
        else:
            n_QP_2D = np.vstack((n_QP_2D, n_QP_1D))

    t = np.arange(0, 41, 1)
    t = t*3*1.21e3
    # n_QP_2D = np.array([d0, d1, d2, d3, d4, d5, d6, d7])
    n_QP_avg = np.average(n_QP_2D, axis=0)
    n_QP_se = np.std(n_QP_2D, axis=0) / np.sqrt(len(n_QP_2D))
    print(repr(n_QP_avg))
    print(repr(n_QP_se))
    # print(n_QP2D)
    plt.xlabel('Poison Length Before T1 Exp (phase slips)')
    plt.ylabel('$n_{qp}$')
    plt.errorbar(t, n_QP_avg, yerr=n_QP_se)

    n = 5  # use the linear part to fit
    [rate, offset] = np.polyfit(t[:n], n_QP_avg[:n], 1)
    plt.plot(t[:n], t[:n] * rate + offset,
             label='QPs/slip={0:.4g}'.format(rate), c='k')
    plt.legend(frameon=False, loc=2, prop={'size': 14})
    plt.show()

if 0:  # QP processed data for single paper plot
    n_QP_avg = np.array([
        0.07380851, 0.1883708, 0.20277472, 0.21447975, 0.30098651,
        0.30824001, 0.30831531, 0.32838444, 0.35551269, 0.39528372,
        0.34032246, 0.35153184, 0.33951708, 0.36688728, 0.39757237,
        0.32735222, 0.37960038, 0.4004139, 0.41697842, 0.4049541,
        0.38849663, 0.32200013, 0.33726678, 0.44671311, 0.39102447,
        0.35429655, 0.39054762, 0.38914731, 0.34114217, 0.4246627,
        0.3550035, 0.36221882, 0.33660496, 0.37218557, 0.39228922,
        0.3995606, 0.32713506, 0.43436966, 0.36348592, 0.36509012,
        0.3987119])
    n_QP_se = np.array([
        0.01963403, 0.02365465, 0.03646389, 0.035712, 0.01705155,
        0.04693009, 0.03483148, 0.03106157, 0.03394227, 0.03039478,
        0.02839453, 0.04214491, 0.03395795, 0.04577062, 0.05306064,
        0.04080479, 0.05910785, 0.03825893, 0.0390526, 0.04640127,
        0.05420951, 0.03128001, 0.0207645, 0.04505482, 0.04009673,
        0.0320043, 0.03328958, 0.04879816, 0.02480708, 0.03648037,
        0.03684882, 0.03801308, 0.03028644, 0.02998843, 0.02728779,
        0.03990929, 0.02155169, 0.04324755, 0.03817371, 0.03885415,
        0.04022622])
    n = 21
    norm_constant = 4.0*1e6
    x_QP_avg = n_QP_avg/norm_constant
    x_QP_se = n_QP_se/norm_constant
    t = np.arange(0, 41, 1)*1.0
    fOver2 = 1.22643
    t_cliff = 2.333*39  # one Clifford gate sequence
    t_cliff_eq = t_cliff*1.22643/1.21/1000  # equivalent Clifford gate sequence us

    # if 1: # just trapping, ignore recombination
    #         # working in n_QP, not x_QP, g=> G=g*n_cp
    t = t*1e-6  # convert time to units of second
    def n_qp(t, G, s):
        return (G/s)*(1-np.exp(-t*s))

    params, covariance = curve_fit(n_qp, t[:n], n_QP_avg[:n])
    std = np.sqrt(np.diag(covariance))  # one standard deviation errors
    [G, s] = params
    n_qp_perPhaseSlip = G/(1.21*3*1e9)
    n_qp_perPhaseSlip_std = std[0]/(1.21*3*1e9)
    s_rate = 1/s
    s_rate_std = (std[1]/s)*s_rate
    print('n_QP per phase slip=', n_qp_perPhaseSlip)
    print('n_QP per phase slip std=', n_qp_perPhaseSlip_std)
    print('s_rate=', s_rate)
    print('s_rate_std=', s_rate_std)

    label_font = 20
    tick_font = 20
    legend_font = 16

    fig = plt.figure(figsize=(8, 5))
    mpl.rc('font', family='Arial')

    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    t_fit = np.arange(0, 20.1, 0.1)
    ax1.plot(t_fit, n_qp(t_fit/1e6, *params), c='k')
    ax1.errorbar(t[:n]*1e6, n_QP_avg[:n], yerr=n_QP_se[:n], fmt="o", color='k', capsize=3.0, ms=5.0)


    ax1.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    ax1.set_xlabel(r'Poisoning Pulse Length ($\mu$s)', fontsize=label_font)
    ax1.set_ylabel('$n_{qp}$', fontsize=label_font)

    ax2.set_xlim(ax1.get_xlim())
    ax2array = np.arange(0, 250, 50)
    ax2.set_xticks(ax2array*t_cliff_eq)
    ax2.set_xticklabels(ax2array)
    ax2.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    ax2.set_xlabel(r"Equivalent Clifford Sequence Length", fontsize=label_font)

    path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
    # plt.savefig(path + '\QPTrapping.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

if 0:  # QP processed data for paper plot with T1 fitting
    """T1 fit starts"""
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = 'T1PoisonSweep'
    experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep_Over4'
    n_QP_2D = np.array([])
    # for i in [0, 1, 2, 3, 4, 5, 6, 7]:
    file_Number = [2]
    T1_2D_file = [file_path.format(date, experiment_name_T1,
                                   experiment_name_T1) + '_{:03d}.mat'.format(
        i) for i in file_Number]
    T1_QP_data = T1_QP_Paper()
    T1_QP_data.add_data_from_matlab(T1_2D_file)
    # n_QP_1D = T1_QP_data.params_2D[:, 0]
    T1_1D_t = T1_QP_data.T1_1D_t/1000
    T1_1D_noPoison = T1_QP_data.T1_1D_noPoison
    T1_1D_noPoisonFit = T1_QP_data.T1_1D_noPoisonFit
    T1_1D_Poison = T1_QP_data.T1_1D_Poison
    T1_1D_PoisonFit = T1_QP_data.T1_1D_PoisonFit

    """T1 fit ends"""

    n_QP_avg = np.array([
        0.07380851, 0.1883708, 0.20277472, 0.21447975, 0.30098651,
        0.30824001, 0.30831531, 0.32838444, 0.35551269, 0.39528372,
        0.34032246, 0.35153184, 0.33951708, 0.36688728, 0.39757237,
        0.32735222, 0.37960038, 0.4004139, 0.41697842, 0.4049541,
        0.38849663, 0.32200013, 0.33726678, 0.44671311, 0.39102447,
        0.35429655, 0.39054762, 0.38914731, 0.34114217, 0.4246627,
        0.3550035, 0.36221882, 0.33660496, 0.37218557, 0.39228922,
        0.3995606, 0.32713506, 0.43436966, 0.36348592, 0.36509012,
        0.3987119])
    n_QP_se = np.array([
        0.01963403, 0.02365465, 0.03646389, 0.035712, 0.01705155,
        0.04693009, 0.03483148, 0.03106157, 0.03394227, 0.03039478,
        0.02839453, 0.04214491, 0.03395795, 0.04577062, 0.05306064,
        0.04080479, 0.05910785, 0.03825893, 0.0390526, 0.04640127,
        0.05420951, 0.03128001, 0.0207645, 0.04505482, 0.04009673,
        0.0320043, 0.03328958, 0.04879816, 0.02480708, 0.03648037,
        0.03684882, 0.03801308, 0.03028644, 0.02998843, 0.02728779,
        0.03990929, 0.02155169, 0.04324755, 0.03817371, 0.03885415,
        0.04022622])
    n = 21
    norm_constant = 4.0*1e6
    x_QP_avg = n_QP_avg/norm_constant
    x_QP_se = n_QP_se/norm_constant
    t = np.arange(0, 41, 1)*1.0
    fOver2 = 1.22643
    t_cliff = 2.333*39  # one Clifford gate sequence
    t_cliff_eq = t_cliff*1.22643/1.21/1000  # equivalent Clifford gate sequence us

    t = t*1e-6  # convert time to units of second
    def n_qp(t, G, s):
        return (G/s)*(1-np.exp(-t*s))

    params, covariance = curve_fit(n_qp, t[:n], n_QP_avg[:n])
    std = np.sqrt(np.diag(covariance))  # one standard deviation errors
    [G, s] = params
    n_qp_perPhaseSlip = G/(1.21*3*1e9)
    n_qp_perPhaseSlip_std = std[0]/(1.21*3*1e9)
    s_rate = 1/s
    s_rate_std = (std[1]/s)*s_rate
    # print('n_QP per phase slip=', n_qp_perPhaseSlip)
    # print('n_QP per phase slip std=', n_qp_perPhaseSlip_std)
    # print('s_rate=', s_rate)
    # print('s_rate_std=', s_rate_std)

    label_font = 20
    tick_font = 20
    legend_font = 18

    fig, axs = plt.subplots(2, figsize=(8, 12)
                            # )
                            ,gridspec_kw = {'hspace': 0.4})
    mpl.rc('font', family='Arial')

    axs[0].scatter(T1_1D_t, T1_1D_noPoison, c='b', label='no poisoning')
    axs[0].plot(T1_1D_t, T1_1D_noPoisonFit, c='b')
    axs[0].scatter(T1_1D_t, T1_1D_Poison, c='r', label='7 $\mu$s poisoning')
    axs[0].plot(T1_1D_t, T1_1D_PoisonFit, c='r')
    axs[0].tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    axs[0].legend(loc=(0.36, 0.56), frameon=False, prop={'size': legend_font}, handletextpad=0)
    axs[0].set_xlabel(r'Delay ($\mu$s)', fontsize=label_font)
    axs[0].set_ylabel('$P_{1}$', fontsize=label_font)
    axs[0].set_xlim([-3, 101])
    axs[0].set_ylim([-0.05, 0.95])
    # axs[0].set_yscale('log')

    ax2 = axs[1].twiny()


    t_fit = np.arange(0, 20.1, 0.1)
    axs[1].plot(t_fit, n_qp(t_fit/1e6, *params), c='k')
    axs[1].errorbar(t[:n]*1e6, n_QP_avg[:n], yerr=n_QP_se[:n], fmt="o", color='k', capsize=3.0, ms=5.0)


    axs[1].tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    axs[1].set_xlabel(r'Poisoning Pulse Length ($\mu$s)', fontsize=label_font)
    axs[1].set_ylabel('$n_{qp}$', fontsize=label_font)
    axs[1].set_xlim([-0.5, 20.5])
    axs[1].set_ylim([-0.02, 0.48])

    ax2.set_xlim(axs[1].get_xlim())
    ax2array = np.arange(0, 250, 50)
    ax2.set_xticks(ax2array*t_cliff_eq)
    ax2.set_xticklabels(ax2array)
    ax2.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    ax2.set_xlabel(r"Equivalent Clifford Sequence Length", fontsize=label_font)

    path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
    plt.savefig(path + '\QPTrapping.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    # plt.savefig(path + '\QPTrapping.pdf', format='pdf', dpi=1200)
    plt.show()

if 0:  # QP processed data
    n_QP_avg = np.array([
        0.07380851, 0.1883708, 0.20277472, 0.21447975, 0.30098651,
        0.30824001, 0.30831531, 0.32838444, 0.35551269, 0.39528372,
        0.34032246, 0.35153184, 0.33951708, 0.36688728, 0.39757237,
        0.32735222, 0.37960038, 0.4004139, 0.41697842, 0.4049541,
        0.38849663, 0.32200013, 0.33726678, 0.44671311, 0.39102447,
        0.35429655, 0.39054762, 0.38914731, 0.34114217, 0.4246627,
        0.3550035, 0.36221882, 0.33660496, 0.37218557, 0.39228922,
        0.3995606, 0.32713506, 0.43436966, 0.36348592, 0.36509012,
        0.3987119])
    n_QP_se = np.array([
        0.01963403, 0.02365465, 0.03646389, 0.035712, 0.01705155,
        0.04693009, 0.03483148, 0.03106157, 0.03394227, 0.03039478,
        0.02839453, 0.04214491, 0.03395795, 0.04577062, 0.05306064,
        0.04080479, 0.05910785, 0.03825893, 0.0390526, 0.04640127,
        0.05420951, 0.03128001, 0.0207645, 0.04505482, 0.04009673,
        0.0320043, 0.03328958, 0.04879816, 0.02480708, 0.03648037,
        0.03684882, 0.03801308, 0.03028644, 0.02998843, 0.02728779,
        0.03990929, 0.02155169, 0.04324755, 0.03817371, 0.03885415,
        0.04022622])
    # print(n_QP2D)
    # normalized to x_{qp}
    norm_constant = 4.0*1e6
    x_QP_avg = n_QP_avg/norm_constant
    x_QP_se = n_QP_se/norm_constant
    t = np.arange(0, 41, 1)*1.0
    # plt.xlabel('Poison Length Before T1 Exp (phase slips)')
    plt.xlabel('SFQ Length ($\mu$s)')
    plt.ylabel('$n_{qp}$')
    # plt.errorbar(t, n_QP_avg, yerr=n_QP_se)
    # plt.errorbar(t, x_QP_avg, yerr=x_QP_se)

    def x_qp(t, x0, gr, r_p, offset):
        return x0*(1-(1-r_p)/(np.exp(t/np.sqrt(4*gr))-r_p))+offset

    if 0:
        def n_qp(t, n0, tau, r_p, offset):
            return n0*(1-(1-r_p)/(np.exp(t/tau)-r_p))+offset

        p0 = [0.3, 2.0, 0.05, 0.01]  # initial guess
        bounds = [(0.1, 1.0, 0, 0), (0.5, 20, 0.2, 0.1)]  # bounds
        params, covariance = curve_fit(n_qp, t, n_QP_avg, p0=p0,
                                           bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        print('params=', params)
        print('std=', std)
        plt.plot(t, n_qp(t, *params))
        plt.errorbar(t, n_QP_avg, yerr=n_QP_se)
        plt.show()

    if 0:
        # def x_qp_full(t, r, s, sqrt_of, r_p, offset):
        #     return ((sqrt_of-s)/(2*r))*(1.0-(1-r_p)/(np.exp(t*sqrt_of)-r_p))+offset
        # def n_qp_full_0(t, n0, sqrt_of, r_p, offset):
        #     return n0*(1.0-(1-r_p)/(np.exp(t*sqrt_of)-r_p))+offset
        def n_qp(t, r, s):
            return ((1.0/np.sqrt(3.897)-s)/(2*r))*(1-0.9/(np.exp(t/3.897)-0.1))+0.085
        # print('x_QP_avg=', x_QP_avg)
        # p0 = [3.897, 0.1 0.085]  # initial guess
        # bounds = [(0.1, 1.0, 0, 0), (0.5, 20, 0.2, 0.1)]  # bounds
        params, covariance = curve_fit(n_qp, t, n_QP_avg)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        print('params=', params)
        print('std=', std)
        plt.plot(t, n_qp(t, *params))
        plt.errorbar(t, n_QP_avg, yerr=n_QP_se)
        # plt.ylim([0, 2e-7])
        plt.show()
    # print()

    if 1: # just trapping, ignore recombination
            # working in n_QP, not x_QP, g=> G=g*n_cp
        t = t*1e-6
        def n_qp(t, G, s):
            return (G/s)*(1-np.exp(-t*s))

        params, covariance = curve_fit(n_qp, t, n_QP_avg)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        [G, s] = params
        print('params=', params)
        print('std=', std)
        n_qp_perPhaseSlip = G/(1.21*3*1e9)
        print('n_QP per phase slip=', n_qp_perPhaseSlip)
        plt.plot(t, n_qp(t, *params))
        plt.errorbar(t, n_QP_avg, yerr=n_QP_se)
        # plt.ylim([0, 2e-7])
        plt.show()


    # p0 = [1e-7, 2.0, 0.05, 0.01*1e-7]  # initial guess
    # bounds = [(1e-8, 1.0, 0, 0), (1e-6, 6.0, 0.2, 1*1e-7)]  # bounds
    # params, covariance = curve_fit(x_qp, t, x_QP_avg, p0=p0,
    #                                    bounds=bounds)
    # std = np.sqrt(np.diag(covariance))  # one standard deviation errors

    if 0: # the linear fit
        n = 5  # use the linear part to fit
        [rate, offset] = np.polyfit(t[:n], n_QP_avg[:n], 1)
        plt.plot(t[:n], t[:n] * rate + offset,
                 label='QPs/slip={0:.4g}'.format(rate), c='k')
        plt.legend(frameon=False, loc=2, prop={'size': 14})
        plt.show()

if 0:  # QP data linear
    #        p0 = [0.2, 6e3, 20e3, 0.95, 0.1]  # initial guess
    #     bounds = [(0, 5e3, 15e3, 0.5, 0), (2, 9e3, 30e3, 1.05, 0.3)]  # bounds
    file_path = (
        'Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
    date = 'T1PoisonSweep'
    experiment_name_T1 = 'T1_SFQ_Poison_Time_Sweep_Over4_LinearPart'
    n_QP_2D = np.array([])
    for i in [0, 1, 2]:
        file_Number = [i]
        T1_2D_file = [file_path.format(date, experiment_name_T1,
                                       experiment_name_T1) + '_{:03d}.mat'.format(
            i) for i in file_Number]
        T1_2D_data = T1_QP_2D()
        T1_2D_data.add_data_from_matlab(T1_2D_file)
        n_QP_1D = T1_2D_data.params_2D[:, 0]
        if len(n_QP_2D) == 0:
            n_QP_2D = n_QP_1D
        else:
            n_QP_2D = np.vstack((n_QP_2D, n_QP_1D))

    t = np.arange(0, 5, 0.25)
    t = t * 3 * 1.21e3
    # n_QP_2D = np.array([d0, d1, d2, d3, d4, d5, d6, d7])
    n_QP_avg = np.average(n_QP_2D, axis=0)
    n_QP_se = np.std(n_QP_2D, axis=0) / np.sqrt(len(n_QP_2D))
    # print(n_QP2D)
    plt.xlabel('Poison Length Before T1 Exp (phase slips)')
    plt.ylabel('$n_{qp}$')
    plt.errorbar(t, n_QP_avg, yerr=n_QP_se)

    n = 20  # use the linear part to fit
    [rate, offset] = np.polyfit(t[:n], n_QP_avg[:n], 1)
    plt.plot(t[:n], t[:n] * rate + offset,
             label='QPs/slip={0:.4g}'.format(rate), c='k')
    plt.legend(frameon=False, loc=2, prop={'size': 14})
    plt.show()

if 0:   # 2D data for basic SFQ-based qubit control
    Ibias_path = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13" \
                "/2022May31SFQControl/Rabi_SFQ_Ibias/MATLABData/Rabi_SFQ_Ibias_000.mat"
    Ibias, IbiasPulseDuration, IbiasOccupation = add_2Ddata_from_matlab(
        Ibias_path, 'SFQ1_I_Bias', 'SFQ_Pulse_Duration', 'Weighted_Occupation')
    Ibias_l = 3
    Ibias_r = 28
    Ibias = -780*Ibias[::-1]
    IbiasOccupation = np.flip(IbiasOccupation, axis=1)
    print(len(Ibias))
    print(len(IbiasOccupation[0]))
    Ibias = Ibias[Ibias_l:Ibias_r]
    IbiasOccupation_new = np.zeros((len(IbiasOccupation), len(Ibias)))
    for i in range(len(IbiasOccupation)):
        IbiasOccupation_new[i] = IbiasOccupation[i][Ibias_l:Ibias_r]
    print(len(Ibias))
    print(len(IbiasOccupation_new[0]))
    print('Ibias=', Ibias)

    Chevron_path = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13" \
                "/2022May31SFQControl/Rabi_SFQ_Chevron/MATLABData/Rabi_SFQ_Chevron_001.mat"
    Freq, FreqPulseDuration, FreqOccupation = add_2Ddata_from_matlab(
        Chevron_path, 'SFQ_Drive_Frequency', 'SFQ_Pulse_Duration', 'Weighted_Occupation')

    Phase_path = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13" \
                "/2022May31SFQControl/SFQ_Y/MATLABData/SFQ_Y_001.mat"
    Phase, PhasePulseDuration, PhaseOccupation = add_2Ddata_from_matlab(
        Phase_path, 'SFQ_Y_Gate_Phase_Offset', 'SFQ_Pulse_Duration', 'Weighted_Occupation')


    # fig, axs = plt.subplots((1, 3), figsize=(11, 8),
    #                         gridspec_kw={'height_ratios': [2, 3], 'hspace': 0.15})
    # fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    mpl.rc('font', family='Arial')
    label_font = 16
    tick_font = 15
    legend_font = 14
    fig = plt.figure(figsize=(20, 4))
    ax1 = plt.subplot(141)
    ax2 = plt.subplot(142)
    ax3 = plt.subplot(143, projection='polar')
    ax4 = plt.subplot(144)

    im = ax1.imshow(IbiasOccupation_new, cmap='jet', aspect="auto", vmin=0, vmax=1,
                    extent=[Ibias[0], Ibias[-1], IbiasPulseDuration[0], IbiasPulseDuration[-1]])
    ax1.set_xlabel(r'$I_{b}$ ($\mu$A)', fontsize=label_font)
    ax1.set_ylabel('SFQ Pulse Duration (ns)', fontsize=label_font)
    ax1.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    ax2.imshow(FreqOccupation, cmap='jet', aspect="auto", vmin=0, vmax=1,
                    extent=[Freq[0], Freq[-1], FreqPulseDuration[0], FreqPulseDuration[-1]])
    ax2.set_xlabel(r'$\omega_{d}/2\pi$ (GHz)', fontsize=label_font)
    ax2.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)


    ip3 = InsetPosition(ax2, [1.05, 0, 1, 1])
    ax3.set_axes_locator(ip3)
    ax3.contourf(Phase, PhasePulseDuration, PhaseOccupation, 100, cmap='jet', vmin=0, vmax=1)
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax3.grid(False)

    ip4 = InsetPosition(ax2, [2.1, 0, 0.05, 1])
    ax4.set_axes_locator(ip4)
    ax4.tick_params(labelsize=tick_font)


    fig.colorbar(im, cax=ax4, ax=[ax1, ax2, ax3])

    # fig.colorbar(im, cax=ax4)
    path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
    plt.savefig(path + '\SFQControl.pdf', format='pdf', bbox_inches='tight', dpi=900)
    plt.show()

if 0:   # Gaussian shape SFQ pulse`
    def gauss_pulse(tau, t):
        v = np.zeros(len(t))   # voltage as function of time
        for gi in range(len(v)):
            v[gi] = 2.0678*np.exp(-t[gi]**2/(2*tau**2))/np.sqrt(2*np.pi*tau**2)
        return v
    def gauss_pulse_f(tau, f):
        v_f = np.zeros(len(f))   # voltage as function of time
        for fi in range(len(f)):
            v_f[fi] = 2.0678*np.exp(-(2*np.pi*tau*f[fi])**2/2)
        return v_f
    # t = np.arange(-10, 10.01, 0.01) #
    # f = np.arange(-1, 1.01, 0.01) #
    t = np.linspace(-10, 10.0, 501) #
    f = np.linspace(-1.0, 1.0, 501) #

    tau_500 = 0.5   #
    tau_1 = 1   #
    tau_2 = 2   #
    tau_5 = 5   #
    v_500 = gauss_pulse(tau_500, t)
    f_500 = gauss_pulse_f(tau_500, f)

    v_1 = gauss_pulse(tau_1, t)
    f_1 = gauss_pulse_f(tau_1, f)

    v_2 = gauss_pulse(tau_2, t)
    f_2 = gauss_pulse_f(tau_2, f)

    v_5 = gauss_pulse(tau_5, t)
    f_5 = gauss_pulse_f(tau_5, f)

    PulseHalfPicosec = []
    PulseOnePicosec = []
    PulseTwoPicosec = []
    for i in range(len(t)):
        dHalfPicosec = [(t[i]+10.0)*10**(-12), v_500[i]]
        dOnePicosec = [(t[i]+10.0)*10**(-12), v_1[i]]
        dTwoPicosec = [(t[i]+10.0)*10**(-12), v_2[i]]
        PulseHalfPicosec.append(dHalfPicosec)
        PulseOnePicosec.append(dOnePicosec)
        PulseTwoPicosec.append(dTwoPicosec)
    print('here')
    np.savetxt('PulseHalfPicosec.txt', PulseHalfPicosec)
    np.savetxt('PulseOnePicosec.txt', PulseOnePicosec)
    np.savetxt('PulseTwoPicosec.txt', PulseTwoPicosec)


    fig, axs = plt.subplots(1, ncols=2, figsize=(8, 4))
        # , ncols=2,
        #                     gridspec_kw={'hspace': 0.5})
    label_font = 20
    tick_font = 20
    legend_font = 14
    axs[0].plot(t, v_500)
    axs[0].plot(t, v_1)
    axs[0].plot(t, v_2)
    # axs[0].plot(t, v_5)
    axs[0].set_xlabel('time (ps)')
    axs[0].set_ylabel('V (mV)')

    axs[1].plot(f, f_500)
    axs[1].plot(f, f_1)
    axs[1].plot(f, f_2)
    # axs[1].plot(f, f_5)
    axs[1].set_xlabel('Frequency (THz)')
    axs[1].set_ylabel('V (mV)')
    plt.show()

if 1:   # Gaussian shape SFQ pulse`
    def gauss_pulse(tau, t):
        v = np.zeros(len(t))   # voltage as function of time
        for gi in range(len(v)):
            v[gi] = 2.0678*np.exp(-t[gi]**2/(2*tau**2))/np.sqrt(2*np.pi*tau**2)
        return v
    def gauss_pulse_f(tau, f):
        v_f = np.zeros(len(f))   # voltage as function of time
        for fi in range(len(f)):
            v_f[fi] = 2.0678*np.exp(-(2*np.pi*tau*f[fi])**2/2)
        return v_f
    t = np.arange(-10, 10.01, 0.01) #
    f = np.arange(-1, 1.01, 0.01) #

    tau_500 = 0.5   #
    tau_1 = 1   #
    tau_2 = 2   #
    tau_5 = 5   #
    v_500 = gauss_pulse(tau_500, t)
    f_500 = gauss_pulse_f(tau_500, f)

    v_1 = gauss_pulse(tau_1, t)
    f_1 = gauss_pulse_f(tau_1, f)

    v_2 = gauss_pulse(tau_2, t)
    f_2 = gauss_pulse_f(tau_2, f)

    v_5 = gauss_pulse(tau_5, t)
    f_5 = gauss_pulse_f(tau_5, f)

    PulseHalfPicosec = []
    PulseOnePicosec = []
    PulseTwoPicosec = []
    for i in range(len(t)):
        dHalfPicosec = [(t[i]+10.0)*10**(-12), v_500[i]]
        dOnePicosec = [(t[i]+10.0)*10**(-12), v_1[i]]
        dTwoPicosec = [(t[i]+10.0)*10**(-12), v_2[i]]
        PulseHalfPicosec.append(dHalfPicosec)
        PulseOnePicosec.append(dOnePicosec)
        PulseTwoPicosec.append(dTwoPicosec)

    fig, ax1 = plt.subplots(figsize=(6, 4))
    label_font = 15
    tick_font = 15
    legend_font = 14

    l, b, h, w = .18, .65, .2, .2
    ax2 = fig.add_axes([l, b, w, h])
    ax2.plot(t, v_2, 'k')
    ax2.axis('off')

    ax1.plot(f, f_500, 'r', linewidth=3, label='$\sigma=0.5$ ps')
    ax1.plot(f, f_1, 'y', linewidth=3, label='$\sigma=1$ ps')
    ax1.plot(f, f_2, 'b', linewidth=3, label='$\sigma=2$ ps')
    ax1.set_xlabel('Frequency (THz)', fontsize=label_font)
    ax1.set_ylabel('Scaled Amplitude', fontsize=label_font)


    ax1.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
    # ax1.legend(True, frameon=True)
    ax1.legend(frameon = False, loc=1, prop={'size': 14})

    path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
    plt.savefig(path + '\GaussianPulse.pdf', format='pdf', bbox_inches='tight', dpi=1200)

    plt.show()

