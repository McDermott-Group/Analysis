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
        CW_freq = 1.2355 # units in GHz
        ps = numberofJJs * CW_freq # phase slips conversion, 3 is the junction number,
        # ps = 1
        plt.figure()
        plt.errorbar(sweep_variable_list*ps, params_2D[:, 0], yerr=std_2D[:, 0], ls='none')
        plt.scatter(sweep_variable_list*ps, params_2D[:, 0])
        if 0:   # recovery rate fitting
            r = 22
            t = sweep_variable_list[:r]
            n_qp = params_2D[:, 0][:r]
            params, std = self._fitToExp(n_qp, t)
            [T, amp, off] = params
            plt.plot(t, self._f_exp_off(t, *params),
                     label='recover time={0:.4g} (ns)'.format(T))
            plt.xlabel('Recover time (ns)')

        if 1:   # n_qp fitting
            ps = sweep_variable_list*9
            rate = 6e-5
            [m, b] = np.polyfit(ps, params_2D[:, 0], 1)
            # plt.plot(ps, ps*m + b, label='QPs/slip={0:.4g}'.format(m), c='k')
            plt.xlabel('SFQ phase slips', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.ylabel('$n_{\mathrm{QP}}$', fontsize=14)
        # plt.xticks([0, 4000, 8000, 12000, 16000, 20000])
        # plt.ylim([0, 2])
        plt.legend(frameon=False, loc=2, prop={'size':14})
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
        sfq_pulse = np.linspace(0, unit*(x_len-1), x_len)
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
        self.SFQ_freq = None    # units is GHz
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
        phase_slips_rate = SFQ_freq*3*1e9   # 3 junctions phase slips
        plt.figure()
        plt.plot(time, occ_fit, 'o-', label='GammaUp={0:.4g}Hz'.format(GammaUp))
        plt.plot(time, occ, 'k-')
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')

        # n_qp_str = 'n_qp={0:.4g} $\pm$ {1:.4g}'.format(n_qp, n_qp_std)
        # T_qp_str = 'T_qp={0:.4g} $\pm$ {1:.4g} ns'.format(T_qp, T_qp_std)
        # T_1r_str = 'T1_r={0:.4g} $\pm$ {1:.4g} ns'.format(T_1r, T_1r_std)
        SFQ_str = 'SFQ Drive Freq = {0:.4g} GHz'.format(SFQ_freq)
        Up_phase_slip_str = 'Induced P1 per phase slip = {0:.4g}'.format(GammaUp/phase_slips_rate)

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