"""
A python library for SFQ project
Author: Chuan-Hong (Vincent) Liu
Date: 2021Dec28
"""
import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.stats import linregress
from scipy.stats import sem
from sklearn.mixture import GaussianMixture
from scipy.optimize import curve_fit


class T1_QP_1D(object):
    """
    This is for extract 1D T1 value from the matlab data
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.QB_Idle_Gate_Time = []
        self.T1 = None
        self.Gamma_fit_parameters = []  # save each file's fitting parameteres
        self.fit_number_list = []
        self.T1_list = []
        self.params = None
        self.err = None
        self.fit_type = 'Exponential'

    def add_data_from_matlab(self, file_path,
                             data_type1='Weighted_Occupation',
                             # data_type1='Projected_Occupation',
                             data_type2='QB_Idle_Gate_Time',
                             fit_type='QP'):
        f_l = file_path[0]
        occ_2D = np.array([])
        Gamma_fit_parameters = np.array([])
        self.QB_Idle_Gate_Time = np.array(
            noiselib.loadmat(f_l)[data_type2])  # unit is ns
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))  # for future average use
        self.fit_type = fit_type

        if fit_type == 'QP':
            fit_params = self.fitToQPDecay(occ_1D)[1]
        Gamma_fit_parameters = np.hstack((Gamma_fit_parameters, fit_params))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])

            if fit_type == 'QP':
                fit_params = self.fitToQPDecay(occ_1D)[1]

            Gamma_fit_parameters = np.hstack(
                (Gamma_fit_parameters, fit_params))

            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters and units"""

        self.occ_1D_avg = occ_1D
        self.Gamma_fit_parameters = Gamma_fit_parameters

    # def fitLinear(self, occ_1D, pts):
    #     time = self.QB_Idle_Gate_Time[pts]
    #     P1 = occ_1D[pts]
    #     params, covariance = curve_fit(self.f_lin, time, P1, p0=[0.9, 30e4])
    #     self.params = params
    #     return params

    # def fitToExponential(self, occ_1D):
    #     """
    #     Extract the time constant, amp, offset and their uncertainties
    #     :param occ_1D:
    #     :return:
    #     """
    #     time = self.QB_Idle_Gate_Time
    #     P1 = occ_1D
    #     p0 =
    #     bounds = [(0.5, 1e4, 0), (1.05, 1e6, 0.4)]
    #     # params, covariance = curve_fit(self.f_exp_off, time, P1,
    #     #                                bounds=bounds, p0=[0.9, 1e5, 0.1])
    #     # params, covariance = curve_fit(self.f_exp_off, time, P1,
    #     #                                p0=[0.9, 2e4, 0.01])
    #     params, covariance = curve_fit(self.f_exp_off, time, P1, p0=p0, bounds=bounds)
    #     self.params = params
    #     return params

    def fitToQPDecay(self, occ_1D):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = self.QB_Idle_Gate_Time
        P1 = occ_1D
        p0 = [0.1, 20e3, 20e3, 0.95, 0.1]  # initial guess
        bounds = [(0, 1e3, 1e3, 0.5, 0), (10, 1e6, 1e6, 1.05, 0.5)]  # bounds
        params, covariance = curve_fit(self.f_QP, time, P1, p0=p0,
                                       bounds=bounds)
        self.params = params
        self.err = np.sqrt(
            np.diag(covariance))  # one standard deviation errors
        # print('params, err=', params, self.err)
        return params

    # def f_lin(self, t, intercept, Gamma):
    #     return intercept - t * Gamma

    # def f_exp(self, t, amp, Gamma):
    #     return amp * np.exp(-t * Gamma) + 0.05
    #
    # def f_exp_off(self, t, amp, Gamma, off):
    #     return amp * np.exp(-t * Gamma) + off

    def f_QP(self, t, n_qp, T_qp, T_1r, amp, off):
        return amp * np.exp(n_qp * (np.exp(-t / T_qp) - 1) - t / T_1r) + off

    def plot_histogram(self):
        Gamma = self.Gamma_fit_parameters
        T1 = [10 ** 6 / g for g in Gamma]  # in us
        plt.hist(T1, bins=10)
        plt.show()

    def plot(self):
        time = self.QB_Idle_Gate_Time
        occ = self.occ_1D_avg
        plt.figure(1)
        plt.plot(time, occ, 'k-')

        [n_qp, T_qp, T_1r, amp, off] = self.params
        [n_qp_std, T_qp_std, T_1r_std, amp_std, off_std] = self.err
        plt.plot(time, self.f_QP(time, n_qp, T_qp, T_1r, amp, off))

        n_qp_str = 'n_qp={0:.4g} $\pm$ {1:.4g}'.format(n_qp, n_qp_std)
        T_qp_str = 'T_qp={0:.4g} $\pm$ {1:.4g} ns'.format(T_qp, T_qp_std)
        T_1r_str = 'T1_r={0:.4g} $\pm$ {1:.4g} ns'.format(T_1r, T_1r_std)

        plt.xlabel('time (ns)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.title(n_qp_str + '\n' + T_qp_str + ', ' + T_1r_str)
        plt.grid()
        plt.legend()
        plt.show()


class T1_QP_2D(object):
    def __init__(self):
        self.P1_2D = None
        self.params_2D = None
        self.std_2D = None

        self.QB_Idle_Gate_Time = None
        self.sweep_variable_list = None
        self.sweep_variable_name = None

    def add_data_from_matlab(self, file_path,
                             # data_type='Weighted_Occupation',
                             data_P1_2D='Projected_Occupation',
                             data_Idle_Time='QB_Idle_Gate_Time',
                             data_sweep_variable='SFQ_Pulse_Duration'
                             # data_sweep_variable='SFQ_Drive_to_QB'
                             ):
        f = file_path[0]
        data = noiselib.loadmat(f)
        P1_2D = np.array(data[data_P1_2D])
        QB_Idle_Gate_Time = np.array(data[data_Idle_Time])
        sweep_variable_list = np.array(data[data_sweep_variable])

        self.P1_2D = P1_2D

        self.QB_Idle_Gate_Time = QB_Idle_Gate_Time
        self.sweep_variable_list = sweep_variable_list
        self.sweep_variable_name = data_sweep_variable

        self._analyze()
        self._plot()
        # sweep_variable_name =
        # print('P1_2D=', P1_2D)
        # print('QB_Idle_Gate_Time=', QB_Idle_Gate_Time)
        # print('sweep_variable_list=', sweep_variable_list)
        return

    def _analyze(self):
        P1_2D = self.P1_2D
        t = self.QB_Idle_Gate_Time
        params_2D = []
        std_2D = []
        ### 1D fit
        for i in range(len(P1_2D)):
            P1_1D = P1_2D[i]
            params, std = self._fitToQPDecay(P1_1D, t)
            params_2D.append(params)
            std_2D.append(std)

        ### 2D data update
        self.params_2D = np.array(params_2D)
        self.std_2D = np.array(std_2D)
        return

    def _plot(self):
        params_2D = self.params_2D
        std_2D = self.std_2D
        sweep_variable_list = self.sweep_variable_list
        sweep_variable_name = self.sweep_variable_name
        # print('params_2D=', params_2D)
        numberofJJs = 3
        CW_freq = 1.2355  # units in GHz
        ps = numberofJJs * CW_freq  # phase slips conversion, 3 is the junction number,
        # ps = 1
        plt.figure()
        plt.errorbar(sweep_variable_list * ps, params_2D[:, 0],
                     yerr=std_2D[:, 0], ls='none')
        plt.scatter(sweep_variable_list * ps, params_2D[:, 0])
        if 0:  # recovery rate fitting
            r = 22
            t = sweep_variable_list[:r]
            n_qp = params_2D[:, 0][:r]
            params, std = self._fitToExp(n_qp, t)
            [T, amp, off] = params
            plt.plot(t, self._f_exp_off(t, *params),
                     label='recover time={0:.4g} (ns)'.format(T))
            plt.xlabel('Recover time (ns)')

        if 1:  # n_qp fitting
            ps = sweep_variable_list * 9
            rate = 6e-5
            [m, b] = np.polyfit(ps, params_2D[:, 0], 1)
            # plt.plot(ps, ps*m + b, label='QPs/slip={0:.4g}'.format(m), c='k')
            plt.xlabel('SFQ phase slips', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.ylabel('$n_{\mathrm{QP}}$', fontsize=14)
        # plt.xticks([0, 4000, 8000, 12000, 16000, 20000])
        # plt.ylim([0, 2])
        plt.legend(frameon=False, loc=2, prop={'size': 14})
        plt.show()
        return

    def _fitToQPDecay(self, P1_1D, time):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = time
        P1 = P1_1D
        p0 = [0.1, 20e3, 20e3, 0.95, 0.1]  # initial guess
        bounds = [(0, 1e3, 1e3, 0.5, 0), (10, 1e6, 1e6, 1.05, 0.5)]  # bounds
        params, covariance = curve_fit(self._f_QP, time, P1, p0=p0,
                                       bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        return list(params), list(std)

    def _f_QP(self, t, n_qp, T_qp, T_1r, amp, off):
        return amp * np.exp(n_qp * (np.exp(-t / T_qp) - 1) - t / T_1r) + off

    def _fitToExp(self, n_qp, time):
        """
        extract recovery rate
        :param n_qp:
        :param time:
        :return:
        """
        time = time
        n_qp = n_qp
        p0 = [20e3, 2, 0]  # initial guess
        bounds = [(1e3, 0, 0), (1e5, 10, 1)]  # bounds
        params, covariance = curve_fit(self._f_exp_off, time, n_qp, p0=p0,
                                       bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        return list(params), list(std)

    def _f_exp_off(self, t, T, amp, off):
        return amp * np.exp(-t / T) + off


class T2_1D(object):
    """
    This is for extract 1D T1 value from the matlab data
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.QB_Idle_Gate_Time = []
        self.T2 = None
        self.fit_parameters = []  # save each file's fitting parameteres
        self.fit_number_list = []
        self.T2_list = []
        self.params = None

    def add_data_from_matlab(self, file_path,
                             data_type1='Weighted_Occupation',
                             # data_type1='Projected_Occupation',
                             data_type2='QB_Idle_Gate_Time'):
        f_l = file_path[0]
        occ_2D = np.array([])
        fit_parameters = np.array([])
        self.QB_Idle_Gate_Time = np.array(
            noiselib.loadmat(f_l)[data_type2])
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        fit_params = self.fitToRamsey(occ_1D, f_l[-7:-4])[3]
        # print('f_l[-7:-4])=', f_l[-7:-4])
        fit_parameters = np.hstack((fit_parameters, fit_params))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            fit_params = self.fitToRamsey(occ_1D, f[-7:-4])[3]
            # print('f[-7:-4])=', f[-7:-4])
            fit_parameters = np.hstack(
                (fit_parameters, fit_params))

        """Update parameters and units"""

        self.occ_1D_avg = occ_1D
        self.fit_parameters = np.array(fit_parameters)

    def fitToRamsey(self, occ_1D, file_number):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param file_number:
        :param occ_1D:
        :return:
        """
        time = self.QB_Idle_Gate_Time
        P1 = occ_1D
        bounds = [(0.3, 1e3, 1e3, 100, 0, 0), (1.05, 1e7, 1e10, 1000, pi, 0.5)]
        p0 = [1, 10e4, 10e4, 300, 0, 0]
        params, covariance = curve_fit(self.f_Ramsey, time, P1, bounds=bounds,
                                       p0=p0)
        self.occ_1D_avg = occ_1D
        self.params = params
        self._plot_fit(file_number, save=True)

        return params

    def f_Ramsey(self, t, amp, T_white, T_corr, T_Ramsey, phase_offset,
                 amp_offset):
        return amp * np.exp(-t / T_white - (t / T_corr) * (t / T_corr)) * \
               np.sin(2 * pi * (t / T_Ramsey) + phase_offset) + amp_offset

    def _plot_fit(self, file_number, save=False):
        time = self.QB_Idle_Gate_Time
        occ = self.occ_1D_avg
        [amp, T_white, T_corr, T_Ramsey, phase_offset,
         amp_offset] = self.params
        plt.figure()
        plt.plot(time, occ, 'k-')
        plt.plot(time, [amp * np.exp(-t / T_white - (t / T_corr) ** 2) * \
                        np.sin(2 * pi * (
                                t / T_Ramsey) + phase_offset) + amp_offset
                        for t in time], '--', label='fit')
        plt.xlabel('time (ns)')
        plt.ylabel('P1')
        plt.grid()
        # amp_str = 'amp={0:.4g} '.format(amp)
        T_white_str = 'T_whilte={0:.4g} (ns) '.format(T_white)
        T_corr_str = 'T_corr={0:.4g} (ns) '.format(T_corr)
        T_Ramsey_str = 'T_Ramsey={0:.4g} (ns)'.format(T_Ramsey)
        plt.title(T_white_str + '\n' + T_corr_str + T_Ramsey_str)
        plt.legend()
        if save:
            plt.savefig(file_number + '.png')
        else:
            plt.show()
        plt.close()

    def plot_Detuning(self):
        Detuning = 1e3 / (self.fit_parameters)  # units in MHz
        x_len = len(Detuning)
        # sfq_pulse = np.linspace(0, 1000, 51) * 3 * 3
        # sfq_pulse = np.linspace(0, 5000, 500) * 1.2355 * 3
        # sfq_pulse = np.linspace(0, 6500, 13)
        # sfq_pulse = np.linspace(0, 500, 2)
        unit = 1000
        sfq_pulse = np.linspace(0, unit * (x_len - 1), x_len)
        plt.figure()
        plt.plot(sfq_pulse, Detuning)
        # plt.xlabel('phase slips (clock time reference)')
        plt.xlabel('Pulse Time (ns)')
        plt.ylabel('Detuning (MHz)')
        # plt.title('2021Dec25_Drive6.3225GHz')
        tt = '2022May19_1.2355GHz_Q1'
        plt.title(tt)
        # plt.savefig(tt + '.png')
        plt.show()


class QP_Up(object):
    """
    This is for analyze the QP up transition rate
    """

    def __init__(self):
        self.occ_1D = []
        self.occ_1D_fit = []
        self.SFQ_freq = None  # units is GHz
        # self.occ_1D_error_band = []
        self.Time_Dep = []
        self.GammaUp = None
        self.GammaUp_std_err = None
        self.temp = None  # units mK
        self.GammaUp_list = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Projected_Occupation',
                             data_type2='SFQ_Pulse_Duration'):
        f_l = file_path[0]
        occ_2D = np.array([])

        # SFQ units is GHz
        SFQ_freq = noiselib.loadmat_ExptVars(f_l)['SFQ_Drive_Frequency']
        # print('SFQ_freq=', SFQ_freq)

        Time_Dep = np.array(noiselib.loadmat(f_l)[data_type2]) * 1e-9
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        Time_Dep = Time_Dep[:]
        occ_1D = occ_1D[:]
        occ_2D = np.hstack((occ_2D, occ_1D))
        #
        # for f in file_path[1:]:
        #     data = noiselib.loadmat(f)
        #     occ_1D = np.array(data[data_type1])
        #     occ_2D = np.vstack((occ_2D, occ_1D))
        """Update parameters and units"""

        self.Time_Dep = Time_Dep
        self.SFQ_freq = SFQ_freq
        self.occ_1D = occ_1D
        self.fitToGammaUp()

    def fitToGammaUp(self):
        time = self.Time_Dep
        P1 = self.occ_1D
        slope, intercept, r_value, p_value, std_err = linregress(time, P1)
        self.GammaUp = slope  # units in Hz
        self.occ_1D_fit = slope * time + intercept
        self.GammaUp_std_err = std_err
        return None

    def plot(self, name='', save=False):
        time = self.Time_Dep * 10 ** 6
        occ = self.occ_1D
        occ_fit = self.occ_1D_fit
        GammaUp = self.GammaUp
        SFQ_freq = self.SFQ_freq  # converted in to GHz
        phase_slips_rate = SFQ_freq * 3 * 1e9  # 3 junctions phase slips
        plt.figure()
        plt.plot(time, occ_fit, 'o-',
                 label='GammaUp={0:.4g}Hz'.format(GammaUp))
        plt.plot(time, occ, 'k-')
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')

        # n_qp_str = 'n_qp={0:.4g} $\pm$ {1:.4g}'.format(n_qp, n_qp_std)
        # T_qp_str = 'T_qp={0:.4g} $\pm$ {1:.4g} ns'.format(T_qp, T_qp_std)
        # T_1r_str = 'T1_r={0:.4g} $\pm$ {1:.4g} ns'.format(T_1r, T_1r_std)
        SFQ_str = 'SFQ Drive Freq = {0:.4g} GHz'.format(SFQ_freq)
        Up_phase_slip_str = 'Induced P1 per phase slip = {0:.4g}'.format(
            GammaUp / phase_slips_rate)

        plt.title(SFQ_str + '\n' + Up_phase_slip_str)

        plt.grid()
        plt.legend()
        plt.show()
        return


class DeepSubharmonics(object):
    """
    This is for average multi SFQ-based rabi to clearly see the steps from
    discretization of SFQ operation
    """

    def __init__(self):
        self.occ_1D = []
        self.occ_1D_avg = []
        self.Time_Dep = []

    def add_data_from_matlab(self, file_path,
                             # data_type1='Projected_Occupation',
                             data_type1='Weighted_Occupation',
                             data_type2='SFQ_Pulse_Duration'):
        f0 = file_path[0]
        data0 = noiselib.loadmat(f0)
        occ_2D = data0[data_type1]
        Time_Dep = data0[data_type2]
        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters"""
        self.Time_Dep = Time_Dep
        self.occ_1D_avg = np.average(occ_2D, axis=0)

    def plot(self, name='', save=False):
        time = self.Time_Dep
        occ = self.occ_1D_avg
        plt.figure()
        plt.plot(time, occ, 'k-')
        plt.xlabel('time (ns)', fontsize=24)
        plt.ylabel('P1', fontsize=24)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        # SFQ_str = 'SFQ Drive Freq = {0:.4g} GHz'.format(SFQ_freq)
        # Up_phase_slip_str = 'Induced P1 per phase slip = {0:.4g}'.format(GammaUp/phase_slips_rate)

        # plt.title(SFQ_str + '\n' + Up_phase_slip_str)

        plt.grid()
        plt.legend()
        plt.show()
        return


class RB(object):
    """
    Analyze the qubit gate fidelity. Import population data, fit, and
    extract the gate fidelity.
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords
        self.NIP = np.array([])  # no interleaved
        self.IP = np.array([])  # interleaved

        self.NIP_avg = np.array([])  # no interleaved avg
        self.NIP_std = np.array([])  # no interleaved standard deviation
        self.NIP_ste = np.array([])  # no interleaved standard error
        self.IP_avg = np.array([])  # interleaved avg
        self.IP_std = np.array([])  # interleaved standard deviation
        self.IP_ste = np.array([])  # interleaved standard error

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             data_type1='No_Interleaved_Probability',
                             data_type2='Interleaved_Probability'):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        self.IntGate = noiselib.loadmat_ExptVars(f0)['Interleaving_Gate']
        NIP_2D = data0[data_type1]
        IP_2D = data0[data_type2]

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            NIP_1D = np.array(data[data_type1])
            IP_1D = np.array(data[data_type2])
            NIP_2D = np.vstack((NIP_2D, NIP_1D))
            IP_2D = np.vstack((IP_2D, IP_1D))

        self.NIP = NIP_2D
        self.IP = IP_2D

    def data_analysis(self):
        self._convert_data()
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _convert_data(self):
        """
        Make the Number of cliffords ascending order and also the
        probability to fidelity
        :return:
        """
        self.NoC = self.NoC
        m = len(self.NoC)
        l1 = len(self.NIP)
        for i in range(l1):
            self.NIP[i] = np.ones(m) - self.NIP[i]
            self.IP[i] = np.ones(m) - self.IP[i]

        self.NIP_avg = np.average(self.NIP, axis=0)
        self.NIP_std = np.std(self.NIP, axis=0)
        self.NIP_ste = self.NIP_std / np.sqrt(l1)
        self.IP_avg = np.average(self.IP, axis=0)
        self.IP_std = np.std(self.IP, axis=0)
        self.IP_ste = self.IP_std / np.sqrt(l1)

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        NIP_avg = self.NIP_avg
        IP_avg = self.IP_avg
        NoC = self.NoC
        p0 = [0.5, 0.985, 0.5]
        params_NIP, covariance_NIP = curve_fit(self._function, NoC, NIP_avg,
                                               p0=p0)
        params_IP, covariance_IP = curve_fit(self._function, NoC, IP_avg,
                                             p0=p0)
        self.params_NIP = params_NIP
        self.params_IP = params_IP
        self.err_NIP = np.sqrt(np.diag(covariance_NIP))
        self.err_IP = np.sqrt(np.diag(covariance_IP))
        # print('params, err=', params_NIP, self.err_NIP)
        # print('params, err=', params_IP, self.err_IP)

    def _function(self, NoC, A, p, B):
        return A * p ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the (non) interleaved seq error and fidelity, and
        find out the interleaved gate's fidelity
        :return:
        """
        p_NIP = self.params_NIP[1]
        p_IP = self.params_IP[1]
        err_NIP = self.err_NIP[1]
        err_IP = self.err_IP[1]

        r_NI = (1 - p_NIP) / 2  # non-interleaved sequence error
        r_I = (1 - p_IP) / 2  # interleaved sequence error

        r_int = (1 - p_IP / p_NIP) / 2.0  # interleaved gate error
        F = 1 - r_int  # interleaved gate fidelity
        F_std = F * np.sqrt((err_NIP / p_NIP) ** 2 + (err_IP / p_IP) ** 2)
        self.r_NI = r_NI
        self.r_I = r_I
        self.F = F
        self.F_std = F_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)
        IntGate = self.IntGate
        NIP_avg = self.NIP_avg
        NIP_ste = self.NIP_ste
        IP_avg = self.IP_avg
        IP_ste = self.IP_ste
        """fitting results"""
        params_NIP = self.params_NIP
        params_IP = self.params_IP
        r_NI = self.r_NI
        r_I = self.r_I
        F = self.F
        F_std = self.F_std

        """real plot"""
        plt.figure()
        plt.errorbar(NoC, NIP_avg, yerr=NIP_ste, fmt="o", color='k')
        plt.plot(NoC_fit, self._function(NoC_fit, *params_NIP), 'k', label=
        'Interleaved Gate: {}, AvgError = {:.3f} %'.format('None', r_NI * 100))

        plt.errorbar(NoC, IP_avg, yerr=IP_ste, fmt="o", color='r')
        plt.plot(NoC_fit, self._function(NoC_fit, *params_IP), 'r', label=
        'Interleaved Gate: {}, AvgError = {:.3f} %'.format(IntGate, r_I * 100))

        plt.title('{} Gate Fidelity = {:.3f} $\pm$ {:.3f} %'.
                  format(IntGate, F * 100, F_std * 100))

        plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Sequence Fidelity')
        plt.show()


class RB_AllGates(object):
    """
    Analyze the qubit gate fidelity. Import population data, fit, and
    extract the gate fidelity. This is for all gates
    ['X', 'X/2', '-X/2', 'Y', 'Y/2', '-Y/2'] with one ref sequence
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords

        self.Ref_sf_2D = np.array([])  # 2D ref sequence fidelity data
        self.Ref_sf_avg = np.array([])  # ref sequence fidelity avg
        self.Ref_sf_std = np.array(
            [])  # ref sequence fidelity standard deviation
        self.Ref_sf_ste = np.array([])  # ref sequence fidelity standard error

        self.X_sf_2D = np.array([])
        self.Y_sf_2D = np.array([])
        self.XOver2_sf_2D = np.array([])
        self.XOver2Neg_sf_2D = np.array([])
        self.YOver2_sf_2D = np.array([])
        self.YOver2Neg_sf_2D = np.array([])

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             data_type1='Ref_Sequence_Fidelity',
                             data_type2='X_Sequence_Fidelity',
                             data_type3='Y_Sequence_Fidelity',
                             data_type4='XOver2_Sequence_Fidelity',
                             data_type5='XOver2Neg_Sequence_Fidelity',
                             data_type6='YOver2_Sequence_Fidelity',
                             data_type7='YOver2Neg_Sequence_Fidelity',
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        Ref_sf_2D = data0[data_type1]
        X_sf_2D = data0[data_type2]
        Y_sf_2D = data0[data_type3]
        XOver2_sf_2D = data0[data_type3]
        XOver2Neg_sf_2D = data0[data_type3]
        YOver2_sf_2D = data0[data_type3]
        YOver2Neg_sf_2D = data0[data_type3]

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            Ref_sf_1D = np.array(data[data_type1])
            Ref_sf_2D = np.vstack((Ref_sf_2D, Ref_sf_1D))
            X_sf_1D = np.array(data[data_type2])
            X_sf_2D = np.vstack((X_sf_2D, X_sf_1D))
            Y_sf_1D = np.array(data[data_type3])
            Y_sf_2D = np.vstack((Y_sf_2D, Y_sf_1D))
            XOver2_sf_1D = np.array(data[data_type4])
            XOver2_sf_2D = np.vstack((XOver2_sf_2D, XOver2_sf_1D))
            XOver2Neg_sf_1D = np.array(data[data_type5])
            XOver2Neg_sf_2D = np.vstack((XOver2Neg_sf_2D, XOver2Neg_sf_1D))
            YOver2_sf_1D = np.array(data[data_type6])
            YOver2_sf_2D = np.vstack((YOver2_sf_2D, YOver2_sf_1D))
            YOver2Neg_sf_1D = np.array(data[data_type7])
            YOver2Neg_sf_2D = np.vstack((YOver2Neg_sf_2D, YOver2Neg_sf_1D))

        self.Ref_sf_2D = Ref_sf_2D
        self.X_sf_2D = X_sf_2D
        self.Y_sf_2D = Y_sf_2D
        self.XOver2_sf_2D = XOver2_sf_2D
        self.XOver2Neg_sf_2D = XOver2Neg_sf_2D
        self.YOver2_sf_2D = YOver2_sf_2D
        self.YOver2Neg_sf_2D = YOver2Neg_sf_2D

    def data_analysis(self):
        self._convert_data()
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _convert_data(self):
        """
        Make 2D raw data to 1D average
        :return:
        """
        self.NoC = self.NoC
        m = len(self.NoC)
        NoS = len(self.Ref_sf_2D)  # number of seqs

        self.Ref_sf_avg = np.average(self.Ref_sf_2D, axis=0)
        self.Ref_sf_std = np.std(self.Ref_sf_2D, axis=0)
        self.Ref_sf_ste = self.Ref_sf_std / np.sqrt(NoS)

        self.X_sf_avg = np.average(self.X_sf_2D, axis=0)
        self.X_sf_std = np.std(self.X_sf_2D, axis=0)
        self.X_sf_ste = self.X_sf_std / np.sqrt(NoS)

        self.Y_sf_avg = np.average(self.Y_sf_2D, axis=0)
        self.Y_sf_std = np.std(self.Y_sf_2D, axis=0)
        self.Y_sf_ste = self.Y_sf_std / np.sqrt(NoS)

        self.XOver2_sf_avg = np.average(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_std = np.std(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_ste = self.XOver2_sf_std / np.sqrt(NoS)

        self.XOver2Neg_sf_avg = np.average(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_std = np.std(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_ste = self.XOver2Neg_sf_std / np.sqrt(NoS)

        self.YOver2_sf_avg = np.average(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_std = np.std(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_ste = self.YOver2_sf_std / np.sqrt(NoS)

        self.YOver2Neg_sf_avg = np.average(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_std = np.std(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_ste = self.YOver2Neg_sf_std / np.sqrt(NoS)
    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        NoC = self.NoC
        p0 = [0.5, 0.985, 0.5]

        Ref_sf_avg = self.Ref_sf_avg
        params_Ref, covariance_Ref = curve_fit(self._function, NoC, Ref_sf_avg,
                                               p0=p0)
        self.params_Ref = params_Ref
        self.err_Ref = np.sqrt(np.diag(covariance_Ref))

        X_sf_avg = self.X_sf_avg
        params_X, covariance_X = curve_fit(self._function, NoC, X_sf_avg, p0=p0)
        self.params_X = params_X
        self.err_X = np.sqrt(np.diag(covariance_X))

        Y_sf_avg = self.Y_sf_avg
        params_Y, covariance_Y = curve_fit(self._function, NoC, Y_sf_avg, p0=p0)
        self.params_Y = params_Y
        self.err_Y = np.sqrt(np.diag(covariance_Y))

        XOver2_sf_avg = self.XOver2_sf_avg
        params_XOver2, covariance_XOver2 = curve_fit(self._function, NoC, XOver2_sf_avg, p0=p0)
        self.params_XOver2 = params_XOver2
        self.err_XOver2 = np.sqrt(np.diag(covariance_XOver2))

        XOver2Neg_sf_avg = self.XOver2Neg_sf_avg
        params_XOver2Neg, covariance_XOver2Neg = curve_fit(self._function, NoC, XOver2Neg_sf_avg, p0=p0)
        self.params_XOver2Neg = params_XOver2Neg
        self.err_XOver2Neg = np.sqrt(np.diag(covariance_XOver2Neg))

        YOver2_sf_avg = self.YOver2_sf_avg
        params_YOver2, covariance_YOver2 = curve_fit(self._function, NoC, YOver2_sf_avg, p0=p0)
        self.params_YOver2 = params_YOver2
        self.err_YOver2 = np.sqrt(np.diag(covariance_YOver2))

        YOver2Neg_sf_avg = self.YOver2Neg_sf_avg
        params_YOver2Neg, covariance_YOver2Neg = curve_fit(self._function, NoC, YOver2Neg_sf_avg, p0=p0)
        self.params_YOver2Neg = params_YOver2Neg
        self.err_YOver2Neg = np.sqrt(np.diag(covariance_YOver2Neg))

    def _function(self, NoC, A, p, B):
        return A * p ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the (non) interleaved seq error and fidelity, and
        find out the interleaved gate's fidelity
        :return:
        """
        p_Ref = self.params_Ref[1]
        err_Ref = self.err_Ref[1]
        r_Ref = (1 - p_Ref) / 2  # non-interleaved sequence error
        self.r_Ref = r_Ref  # avg gate error

        p_X = self.params_X[1]
        err_X = self.err_X[1]
        r_X = (1 - p_X) / 2  # interleaved sequence error
        F_X = (1 + p_X / p_Ref) / 2.0  # X gate fidelity
        F_X_std = F_X * np.sqrt((err_Ref / p_Ref) ** 2 + (err_X / p_X) ** 2)
        self.r_X = r_X  # avg gate+X error
        self.F_X = F_X  # X gate fidelity
        self.F_X_std = F_X_std  # X gate fidelity std

        p_Y = self.params_Y[1]
        err_Y = self.err_Y[1]
        r_Y = (1 - p_Y) / 2  # interleaved sequence error
        F_Y = (1 + p_Y / p_Ref) / 2.0  # X gate fidelity
        F_Y_std = F_Y * np.sqrt((err_Ref / p_Ref) ** 2 + (err_Y / p_Y) ** 2)
        self.r_Y = r_Y  # avg gate+X error
        self.F_Y = F_Y  # X gate fidelity
        self.F_Y_std = F_Y_std  # X gate fidelity std

        p_XOver2 = self.params_XOver2[1]
        err_XOver2 = self.err_XOver2[1]
        r_XOver2 = (1 - p_XOver2) / 2
        F_XOver2 = (1 + p_XOver2 / p_Ref) / 2.0
        F_XOver2_std = F_XOver2 * np.sqrt((err_Ref / p_Ref) ** 2 + (err_XOver2 / p_XOver2) ** 2)
        self.r_XOver2 = r_XOver2
        self.F_XOver2 = F_XOver2
        self.F_XOver2_std = F_XOver2_std

        p_XOver2Neg = self.params_XOver2Neg[1]
        err_XOver2Neg = self.err_XOver2Neg[1]
        r_XOver2Neg = (1 - p_XOver2Neg) / 2
        F_XOver2Neg = (1 + p_XOver2Neg / p_Ref) / 2.0
        F_XOver2Neg_std = F_XOver2Neg * np.sqrt((err_Ref / p_Ref) ** 2 + (err_XOver2Neg / p_XOver2Neg) ** 2)
        self.r_XOver2Neg = r_XOver2Neg
        self.F_XOver2Neg = F_XOver2Neg
        self.F_XOver2Neg_std = F_XOver2Neg_std

        p_YOver2 = self.params_YOver2[1]
        err_YOver2 = self.err_YOver2[1]
        r_YOver2 = (1 - p_YOver2) / 2
        F_YOver2 = (1 + p_YOver2 / p_Ref) / 2.0
        F_YOver2_std = F_YOver2 * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2 / p_YOver2) ** 2)
        self.r_YOver2 = r_YOver2
        self.F_YOver2 = F_YOver2
        self.F_YOver2_std = F_YOver2_std

        p_YOver2Neg = self.params_YOver2Neg[1]
        err_YOver2Neg = self.err_YOver2Neg[1]
        r_YOver2Neg = (1 - p_YOver2Neg) / 2
        F_YOver2Neg = (1 + p_YOver2Neg / p_Ref) / 2.0
        F_YOver2Neg_std = F_YOver2Neg * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2Neg / p_YOver2Neg) ** 2)
        self.r_YOver2Neg = r_YOver2Neg
        self.F_YOver2Neg = F_YOver2Neg
        self.F_YOver2Neg_std = F_YOver2Neg_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 96, 1)

        """real plot"""
        plt.figure()
        plt.errorbar(NoC, self.Ref_sf_avg, yerr=self.Ref_sf_ste, fmt="o", color='k')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Ref), 'k',
                 label='Ref, AvgError = {:.3f} %'.format(self.r_Ref * 100))

        plt.errorbar(NoC, self.X_sf_avg, yerr=self.X_sf_ste, fmt="o", color='r')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_X), 'r',
                 label='X, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_X * 100, self.F_X * 100, self.F_X_std * 100))

        plt.errorbar(NoC, self.Y_sf_avg, yerr=self.Y_sf_ste, fmt="o", color='b')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Y), 'b',
                 label='Y, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_Y * 100, self.F_Y * 100, self.F_Y_std * 100))

        plt.errorbar(NoC, self.XOver2_sf_avg, yerr=self.XOver2_sf_ste, fmt="o", color='y')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2), 'y',
                 label='X/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_XOver2 * 100, self.F_XOver2 * 100, self.F_XOver2_std * 100))

        plt.errorbar(NoC, self.XOver2Neg_sf_avg, yerr=self.XOver2Neg_sf_ste, fmt="o", color='g')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2Neg), 'g',
                 label='-X/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_XOver2Neg * 100, self.F_XOver2Neg * 100, self.F_XOver2Neg_std * 100))

        plt.errorbar(NoC, self.YOver2_sf_avg, yerr=self.YOver2_sf_ste, fmt="o",
                     color='r')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2), 'r',
                 label='Y/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_YOver2 * 100, self.F_YOver2 * 100,
                         self.F_YOver2_std * 100))

        plt.errorbar(NoC, self.YOver2Neg_sf_avg, yerr=self.YOver2Neg_sf_ste,
                     fmt="o", color='r')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2Neg), 'r',
                 label='-Y/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_YOver2Neg * 100, self.F_YOver2Neg * 100,
                         self.F_YOver2Neg_std * 100))

        plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Sequence Fidelity')
        plt.show()
