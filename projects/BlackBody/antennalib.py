"""
A library for the antenna model
"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.stats import linregress
from scipy.stats import sem
from sklearn.mixture import GaussianMixture

from scipy.optimize import curve_fit

"""First import matlab data to python"""


class QP_Up(object):
    """
    This is for analyze the QP up transition rate with premeasurement as
    initialization.
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.Pre_to_RO_Time = []
        self.GammaUp = None
        self.GammaUp_std_err = None
        self.temp = None  # units mK
        self.Chi_nu_Square_list = []
        self.fit_number_list = []
        self.GammaUp_list = []

    def add_data_from_matlab(self, file_path, temp=None,
                             data_type1='Weighted_Occupation_Filtered',
                             # data_type1='Single_Shot_Occupation_Filtered',
                             data_type2='RO_Pre_Drive_to_RO'):
        f_l = file_path[0]
        weights = []
        residual_reps = noiselib.loadmat_ExptVars(f_l)['Residual_Reps']
        weights.append(residual_reps)
        occ_2D = np.array([])
        self.Pre_to_RO_Time = np.array(
            noiselib.loadmat(f_l)[data_type2]) * 1e-9
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))

        for f in file_path[1:]:
            r_rep = noiselib.loadmat_ExptVars(f)['Residual_Reps']
            residual_reps += r_rep
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))
            weights.append(r_rep)

        # print('weights=', weights)
        # print('occ_2D=', occ_2D)
        # print(occ_2D.mean(axis=0))
        """Update parameters and units"""
        occ_1D_avg, error = self._weighted_avg_and_std_2D(occ_2D, weights)

        self.occ_1D_avg = occ_1D_avg
        self.occ_1D_avg_error_band = error
        self.fitToGammaUp()
        if isinstance(temp, str):
            temp = float(temp[:3])
            # print('temp=', temp)
        self.temp = temp

    def _weighted_avg_and_std_2D(self, occ_2D, weights):
        """
        calculate the weighted average and error
        :param occ_2D: [[P1(t0), ..., P1(tf)], [P2(t0), ..., P2(tf)],..., [Pn(t0), ..., Pn(tf)]]
        :param weights: [w1, w2, ..., wn]
        :return: occ_1D_avg, error
        """
        occ_1D_avg = []
        error = []
        for i in range(len(occ_2D[0])):
            p1 = occ_2D[:,
                 i]  # a list of Poplation at the same time [P1(ti), P2(ti), ...,Pn(ti)]
            avg, std = self._weighted_avg_and_std(p1, weights)
            occ_1D_avg.append(avg)
            error.append(std)
        # to get the error in the mean
        error = error / np.sqrt(len(weights))
        return np.array(occ_1D_avg), np.array(error)

    def _weighted_avg_and_std(self, values, weights):
        """
        Return the weighted average and standard deviation.

        values, weights -- Numpy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        variance = np.average((values - average) ** 2, weights=weights)
        return average, np.sqrt(variance)

    def fitToGammaUp(self):
        time = self.Pre_to_RO_Time
        P1 = self.occ_1D_avg
        if False:  # update length of the fit, use the first 25% of the data to fit
            length = int(len(time) / 4)
            time = time[:length]
            P1 = P1[:length]
        slope, intercept, r_value, p_value, std_err = linregress(time, P1)
        # print('std_err=', std_err)
        self.GammaUp = slope  # units in Hz
        self.occ_1D_avg_fit = slope * time + intercept
        self.GammaUp_std_err = std_err
        return None

    def fitToGammaUp_Chi(self):
        """
        sweep the fitting length and get chi^2 to see if the fit is good or not
        :return:
        """
        time = self.Pre_to_RO_Time
        P1 = self.occ_1D_avg
        P1_sig = self.occ_1D_avg_error_band
        l = len(P1)
        Chi_nu_Square_list = []
        fit_number_list = []
        GammaUp_list = []

        for n in range(3, l):
            Chi_Square = 0
            Chi_nu_Square = 0
            slope, intercept, r_value, p_value, std_err = linregress(time[:n],
                                                                     P1[:n])
            P1_fit = time * slope + intercept
            # print('P1_fit[:n]=', P1_fit[:n])
            for i in range(n):
                Chi_Square = Chi_Square + ((P1_fit[i] - P1[i]) ** 2) / (
                        P1_sig[i] ** 2)
            Chi_nu_Square = Chi_Square / (n - 2)
            Chi_nu_Square_list.append(Chi_nu_Square)
            fit_number_list.append(n)
            GammaUp_list.append(slope)  # units in Hz
        self.Chi_nu_Square_list = Chi_nu_Square_list
        self.fit_number_list = fit_number_list
        self.GammaUp_list = GammaUp_list
        return None

    def plot_Chi(self):
        self.fitToGammaUp_Chi()
        Chi_nu_Square_list = self.Chi_nu_Square_list
        fit_number_list = self.fit_number_list
        GammaUp_list = self.GammaUp_list
        time = self.Pre_to_RO_Time * 10 ** 6
        plt.figure(1)
        # plt.plot(fit_number_list, Chi_nu_Square_list, 'k-', label='points')
        plt.plot(time[3:], Chi_nu_Square_list, 'k-', label='Chi_nu')
        plt.xlabel('time (us)')
        plt.ylabel('Chi_Square/nu')
        plt.grid()
        plt.legend()

        plt.figure(2)
        plt.xlabel('time (us)')
        plt.ylabel('Gamma Up (Hz)')
        plt.grid()
        plt.legend()
        plt.plot(time[3:], GammaUp_list, 'k-', label='GammaUp')
        #plt.show(block=False)

    def plot(self,name='',save=False):
        time = self.Pre_to_RO_Time * 10 ** 6
        occ = self.occ_1D_avg
        occ_fit = self.occ_1D_avg_fit
        error = self.occ_1D_avg_error_band
        GammaUp = self.GammaUp
        temp = str(self.temp) + 'mK'
        if not save:
            plt.figure(3)
        plt.plot(time, occ_fit, 'o-',
                 label=temp + '_GammaUp={0:.4g} Hz'.format(GammaUp))
        plt.plot(time, occ, 'k-', label=temp)
        plt.fill_between(time, occ - error, occ + error, alpha=0.2)
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        if save:
            plt.savefig(name)
            plt.show()
        else:
            plt.show()


class Up_array(object):
    """
    This is for collecting Temp, Up_rate, Up_rate_error array
    and plot them
    """

    def __init__(self):
        # self.temp = np.array([])    # units mK
        self.temp = []  # units mK
        self.GammaUp = []  # units Hz
        self.GammaUp_std_err = []  # units Hz

    def insert_data(self, array_tempUpUpError):
        temp = np.append(self.temp, array_tempUpUpError[0])
        GammaUp = np.append(self.GammaUp, array_tempUpUpError[1])
        GammaUp_std_err = np.append(self.GammaUp_std_err,
                                    array_tempUpUpError[2])

        zipped_lists_temp = zip(temp, temp)
        print('zipped_lists_temp=', zipped_lists_temp)
        sorted_zipped_lists_temp = sorted(zipped_lists_temp)
        print('sorted_zipped_lists_temp=', sorted_zipped_lists_temp)
        sorted_temp = np.array(
            [element for _, element in sorted_zipped_lists_temp])
        print('sorted_temp=', sorted_temp)

        zipped_lists_GammaUp_std_err = zip(temp, GammaUp_std_err)
        sorted_zipped_lists_GammaUp_std_err = sorted(
            zipped_lists_GammaUp_std_err)
        sorted_temp_std_err = np.array(
            [element for _, element in sorted_zipped_lists_GammaUp_std_err])

        zipped_lists_GammaUp = zip(temp, GammaUp)
        sorted_zipped_lists_GammaUp = sorted(zipped_lists_GammaUp)
        sorted_temp = np.array(
            [element for _, element in sorted_zipped_lists_GammaUp])

        self.temp = sorted_temp
        self.GammaUp = sorted_GammaUp
        self.GammaUp_std_err = sorted_GammaUp_std_err

    def plot(self):
        #
        Temp = self.temp * 10 ** (-3)
        Gamma = self.GammaUp
        GammaErrorMinus = self.GammaUp - self.GammaUp_std_err
        GammaErrorPlus = self.GammaUp + self.GammaUp_std_err
        # plt.scatter([1 / T for T in T_Pre],
        #             [np.log(1 / G) for G in Gamma_Pre_Q3], label='Q3_Up_Pre')
        # plt.plot([1/T for T in T_Pre_Q4], [np.log(1/G) for G in Gamma_Pre_Q4], label='Q4_Up_Pre')
        # plt.plot([1 / T for T in T_Pre], [0.2 / T - 8.7 for T in T_Pre],
        #          label='Fit, f0=4GHz')
        plt.plot([1 / T for T in Temp], [np.log(1 / G) for G in Gamma],
                 label='Up_Pre')

        plt.fill_between([1 / T for T in Temp],
                         [np.log(1 / G) for G in GammaErrorMinus],
                         [np.log(1 / G) for G in GammaErrorPlus], alpha=0.2)

        # plt.plot([1 / T for T in Temp], [0.2 / T - 8.7 for T in Temp],
        #          label='Fit, f0=4GHz')

        plt.xlabel('1/Temp(Kelvin)')
        plt.ylabel('ln(df/Gamma)')
        plt.grid()
        plt.legend()
        plt.show()


class T1(object):
    """
    This is for extract T1 value from the matlab data
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.QB_Idle_Gate_Time = []
        self.T1 = None
        self.Gamma_fit_parameters = []  # save each file's fitting parameteres
        self.Gamma_fit_parameters = []  # save each file's fitting parameteres
        self.temp = None  # units mK
        self.fit_number_list = []
        self.T1_list = []

    def add_data_from_matlab(self, file_path, temp,
                             data_type1='Weighted_Occupation',
                             # data_type1='Projected_Occupation',
                             data_type2='QB_Idle_Gate_Time',
                             fit_type='Exponential', num_points=5):
        f_l = file_path[0]
        occ_2D = np.array([])
        Gamma_fit_parameters = np.array([])
        self.QB_Idle_Gate_Time = np.array(
            noiselib.loadmat(f_l)[data_type2]) * 1e-9
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))  # for future average use
        # fit_params = np.array([self.fitToExponential(occ_1D)])
        fit_params = self.fitToExponential(occ_1D)[1]
        Gamma_fit_parameters = np.hstack((Gamma_fit_parameters, fit_params))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])

            if fit_type == 'Exponential':
                fit_params = self.fitToExponential(occ_1D)[1]
            else:  # fit_type='Linear'
                fit_params = self.fitLinear(occ_1D, slice(0, num_points))[1]

            Gamma_fit_parameters = np.hstack(
                (Gamma_fit_parameters, fit_params))

            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters and units"""
        # occ_1D_avg, error = self._weighted_avg_and_std_2D(occ_2D, weights)
        # print('Gamma_fit_parameters=', Gamma_fit_parameters)

        self.occ_1D_avg = occ_1D
        self.temp = temp
        self.Gamma_fit_parameters = Gamma_fit_parameters

    def fitLinear(self, occ_1D, pts):
        time = self.QB_Idle_Gate_Time[pts]
        P1 = occ_1D[pts]
        params, covariance = curve_fit(self.f_lin, time, P1, p0=[0.9, 30e4])
        self.params = params
        return params

    def fitToExponential(self, occ_1D):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = self.QB_Idle_Gate_Time
        P1 = occ_1D
        # params, covariance = curve_fit(self.f_exp, time, P1, p0=[0.9, 30e4])
        bounds = [(0.5, 1e4, 0.008), (1.05, 1e6, 0.4)]
        params, covariance = curve_fit(self.f_exp_off, time, P1,
                                       bounds=bounds, p0=[0.9, 1e5, 0.1])
        # params, covariance = curve_fit(self.f_exp, time, P1, p0=[0.9, 3e4])
        self.params = params
        return params

    def f_lin(self, t, intercept, Gamma):
        return intercept - t * Gamma

    def f_exp(self, t, amp, Gamma):
        return amp * np.exp(-t * Gamma) + 0.05

    def f_exp_off(self, t, amp, Gamma, off):
        return amp * np.exp(-t * Gamma) + off

    def plot_histogram(self):
        Gamma = self.Gamma_fit_parameters
        T1 = [10 ** 6 / g for g in Gamma]  # in us
        plt.hist(T1, bins=10)
        plt.show()

    def plot(self):
        time = self.QB_Idle_Gate_Time * 10 ** 6
        occ = self.occ_1D_avg
        amp, Gamma, off = self.params[0], self.params[1], self.params[2]
        temp = str(self.temp) + 'mK'
        plt.figure(1)
        plt.plot(time, occ, 'k-', label=temp)
        plt.plot(time, [amp * np.exp(-t * Gamma / (10 ** 6)) for t in time],
                 '--', label='fit')
        # plt.plot(time, [amp*np.exp(-t*Gamma/(10**6))+off for t in time], '--', label='fit')
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        plt.yscale('log')
        plt.grid()
        amp_str = 'amp={0:.4g} '.format(amp)
        off_str = 'offset={0:.4g}'.format(off)
        T1_str = 'T1={0:.4g} us'.format(1e6 / Gamma)
        plt.title(T1_str + '\n' + amp_str + off_str)
        plt.legend()
        plt.show()


class P1(object):
    """
    This is for extract P1 value from the matlab data. Collaborate with T1
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_std = None
        self.temp = None  # units mK

    def add_data_from_matlab(self, file_path, temp,
                             # data_type1='Weighted_Occupation'):
                             data_type1='Projected_Occupation',
                             data_type2='Projected_Occupation_Std'):
        occ_1D_avg = []
        occ_1D_std = []
        for f in file_path:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_std = np.array(data[data_type2])
            occ_1D_avg.append(np.mean(occ_1D))
            # for i in range(len(occ_1D)):
            #     occ_1D_avg.append(occ_1D[i])
            #     occ_1D_std.append(occ_std[i])

        """Update parameters"""
        # print('occ_1D_avg=', occ_1D_avg)

        self.occ_1D_avg = np.array(occ_1D_avg)
        # self.occ_std = np.mean(occ_1D_std)
        self.temp = temp


class AntennaCoupling(object):
    """
    This is analyze the CST simulation
    """

    def __init__(self):
        self.Antenna = {
            "f": [],
            "Z_Re": [],
            "Z_Im": [],
            "Z_rad": [],
            "Gamma": [],
        }

        self.Junction = {
            "R": None,
            "L": None,
            "C": None,
            "Area": None,
            "Z_j": None,
            "Ic": None
        }
        Gamma = None
        e_c = None
        Ic_f = None   # as a function of freq
        e_c_dB = None
        p_g = None  # photon generation rate
        ref = None  # Al sample box reflection array
        C_eff = 150*1e-21    # F/nm^2

    def import_data(self, file, JJ, C_eff=150*1e-21):
        self.C_eff = C_eff
        self._add_data_from_txt(file)
        self._JJ_Update(JJ)
        self._get_e_c()
        self._get_Ic()
        self._get_p_g()
        self._get_ref()

    def _add_data_from_txt(self, file):
        f = np.loadtxt(file, usecols=[0], skiprows=3)
        Z_Re = np.loadtxt(file, usecols=[1], skiprows=3)
        Z_Im = np.loadtxt(file, usecols=[2], skiprows=3)
        Z_rad = Z_Re + 1j * Z_Im
        self.Antenna["f"] = f * 1e9
        self.Antenna["Z_Re"] = Z_Re
        self.Antenna["Z_Im"] = Z_Im
        self.Antenna["Z_rad"] = Z_rad

    # def _JJ_Update(self, JJ):
    #     R = JJ[0]
    #     L = JJ[1]
    #     C = JJ[2]
    #     A = JJ[3]  # nm*nm
    #     Ic = (np.pi/4)*(380*(1e-6)/R)
    #     C_eff = self.C_eff
    #     omega = 2 * pi * self.Antenna["f"]
    #     if not C:
    #         # print('here!')
    #         C = A * C_eff
    #         tau = R * C
    #         # print('tau=', tau)
    #         Z_j = [(1 - 1j * w * tau) / (1 + w ** 2 * tau ** 2) * R for w in
    #                omega]
    #     else:
    #         print('Q3')
    #         Z_j = [-1j / (w * C - 1 / (w * L)) for w in omega]
    #     self.Junction["R"] = R
    #     self.Junction["L"] = L
    #     self.Junction["C"] = C
    #     self.Junction["A"] = A
    #     self.Junction["Z_j"] = Z_j
    #     self.Junction["Ic"] = Ic

    def _JJ_Update(self, JJ):   # Wirebond included
        R = JJ[0]
        L = JJ[1]
        C = JJ[2]
        A = JJ[3]  # nm*nm
        Ic = (np.pi/4)*(380*(1e-6)/R)
        C_eff = self.C_eff
        omega = 2 * pi * self.Antenna["f"]
        # print('here!')
        C = A * C_eff
        tau = R * C
        # print('tau=', tau)
        length_wb_bias = 0.003   # units mm
        L_wb_bias = mu_0*length_wb_bias   #wb = wirebond
        C_wb_bias = epsilon_0*length_wb_bias    #wb = wirebond
        length_wb_gnd = 0.0003   # units mm
        L_wb_gnd = mu_0*length_wb_gnd   #wb = wirebond
        C_wb_gnd = epsilon_0*length_wb_gnd    #wb = wirebond
        Z_j = []
        for w in omega:
            Z_JJ = 1/(1/R+1j*w*C)
            # Z_wirebonda_gnd = 1/(1/1j*w*L_wb_gnd+1j*w*C_wb_gnd)
            # Z_wirebonda_bias = 1/(1/1j*w*L_wb_bias+1j*w*C_wb_bias)
            # Z_wirebond = 0
            # Z = 1/(1/(Z_JJ + Z_wirebonda_gnd)+1/(Z_wirebonda_bias+50))
            Z = 1/(1/R+1j*w*C)    # no wirebond included
            Z_j.append(Z)
        self.Junction["R"] = R
        self.Junction["L"] = L
        self.Junction["C"] = C
        self.Junction["A"] = A
        self.Junction["Z_j"] = Z_j
        self.Junction["Ic"] = Ic

    def _get_e_c(self):
        Z_j = self.Junction["Z_j"]
        Z_rad = self.Antenna["Z_rad"]
        Gamma = []
        for i in range(len(Z_rad)):
            Gamma.append((Z_rad[i] - np.conj(Z_j[i])) / (Z_rad[i] + Z_j[i]))
        e_c = 1 - (np.abs(Gamma)) ** 2
        e_c_dB = 10 * np.log10(e_c)
        self.Gamma = Gamma
        self.e_c = e_c
        self.e_c_dB = e_c_dB

    # def _get_p_g(self):
    #     e_c = self.e_c
    #     f = self.Antenna["f"]
    #     Ic = self.Junction["Ic"]
    #     R = self.Junction["R"]
    #     p_g = []
    #
    #     P = Ic**2*R
    #     Pf = P * e_c
    #     for i in range(len(f)):
    #         p_g_f = Pf[i]/(h*f[i])
    #         p_g.append(p_g_f)
    #
    #     # print("p_g_f=", p_g_f)
    #     # print("P=", P)
    #     # print("p_g[100:110]", p_g[100:110])
    #     # print("f[100:110]", f[100:110])
    #     # print(len(p_g))
    #     self.p_g = p_g

    def _get_Ic(self):
        """
        get Ic as a function of radiation freq, consider critical current
        as a function of temperature temp is a function of power dissipated on chip
        :return:
        """
        Omega = 40*1e-6 * 1.6*1e-6 * 50*1e-15   # volume of the junction
        Sigma = 2.4*1e9 # the single parameter from Wellstood and Clarke's 1994 paper

        Ic_f = []
        T_f = []
        # Vb = []
        P_heat_f = []

        f = self.Antenna["f"]
        R = self.Junction["R"]*1.0
        phi_0 = h/(2*e)
        # print('R=', R)

        # for fi in f:
        #     vb = fi * phi_0  # convert photon frequency to voltage bias
        #     P_heat = vb**2.0 / R
        #     T = (P_heat/(Omega*Sigma))**(0.2)   # cannot use 1/5, bust use 0.2
        #     T_f.append(T)
        #     ic = self._getIcFromTemp(T)
        #     # Vb.append(vb)
        #     P_heat_f.append(P_heat)
        #     Ic_f.append(ic)

        for fi in f:    # Tinkham
            vb = fi * phi_0  # convert photon frequency to voltage bias
            T_b = 10e-3 # base temp 20mK
            if vb <= 0:
                vb = 0
            T = np.sqrt((T_b**2+3*(e*vb/(2*pi*k))))
            T_f.append(T)
            ic = self._getIcFromTemp(T)
            # Vb.append(vb)
            Ic_f.append(ic)

        # print("Vb=", Vb)
        # print("T_f=", T_f[0:10])
        # print("T_f=", T_f[200:210])
        # print("T_f=", T_f[900:910])
        # print("Ic_f=", Ic_f[0:10])
        # print("Ic_f=", Ic_f[200:210])
        # print("Ic_f=", Ic_f[900:910])
        # print("P_heat_f=", P_heat_f)

        self.Ic_f = Ic_f

    # def _get_p_g(self):
    #     """
    #     photons generated,
    #     :return:
    #     """
    #     e_c = self.e_c
    #     f = self.Antenna["f"]
    #     Ic = self.Junction["Ic"]
    #     R = self.Junction["R"]
    #     p_g = []
    #
    #     P = Ic**2*R
    #     Pf = P * e_c
    #     for i in range(len(f)):
    #         p_g_f = Pf[i]/(h*f[i])
    #         p_g.append(p_g_f)
    #
    #     # print("p_g_f=", p_g_f)
    #     # print("P=", P)
    #     # print("p_g[100:110]", p_g[100:110])
    #     # print("f[100:110]", f[100:110])
    #     # print(len(p_g))
    #     self.p_g = p_g

    def _get_p_g(self):
        """
        photons generated,
        :return:
        """
        e_c = self.e_c
        f = self.Antenna["f"]
        Ic = self.Junction["Ic"]
        Ic_f = self.Ic_f
        R = self.Junction["R"]
        p_g = []

        # P = Ic**2*R
        # Pf = P * e_c
        for i in range(len(f)):
            P = Ic_f[i]**2 * e_c[i] * R
            p_g_f = P/(h*f[i])
            p_g.append(p_g_f)

        self.p_g = p_g

    def _get_ref(self):
        """
        reflection coefficient
        :return:
        """
        f = self.Antenna["f"]
        omega = 2*np.pi*f
        sigma = 1e9
        Z0 = 377    # vacuum impedance
        # mu_0
        # print("mu_0=", mu_0)
        Gamma = []
        Gamma_square = []
        for i in range(len(omega)):
            Z_Al = (1+1j)*np.sqrt((omega[i]*mu_0)/(2*sigma))
            gamma = (Z_Al-Z0)/(Z_Al+Z0)
            Gamma.append(gamma)
            Gamma_square.append(np.abs(gamma))

        # print(Gamma_square)
        self.ref = Gamma_square

        # Gamma = []
        # for i in range(len(Z_rad)):
        #     Gamma.append((Z_rad[i] - np.conj(Z_j[i])) / (Z_rad[i] + Z_j[i]))

    def _getIcFromTemp(self, T):
        """
        Get critical current from temperature
        :param T:
        :return:
        """
        Rn = self.Junction["R"]
        Delta = 190.0*1e-6*e
        IcRn = (pi/4)*(2*Delta/e)*np.tanh(Delta / (2 * k * T))
        Ic = IcRn/Rn
        # if T > 0.202 and T < 0.203:
        #     print('T=', T)
        #     print('IcRn=', IcRn)
        #     print('Delta / (2 * k * T)=', Delta / (2 * k * T))
        #     print('Ic=', Ic)
        return Ic

    def plot(self):
        f = self.Antenna["f"]
        e_c_dB = self.e_c_dB
        plt.plot(f, e_c_dB)
        plt.show()

    # def plot_Q3Mode(self):
    #     f = self.Antenna["f"]
    #     Z_j = self.Junction["Z_j"]
    #     Z_rad = self.Antenna["Z_rad"]
    #     Y_j = [1 / Z for Z in Z_j]
    #     Y_rad = [1 / Z for Z in Z_rad]
    #     Y_tot = np.add(Y_j, Y_rad)
    #     plt.plot(f, Y_tot.real, label='Re')
    #     plt.plot(f, Y_tot.imag, label='Im')
    #     plt.grid()
    #     plt.legend()
    #     plt.show()


class P1_JSweep(object):
    """
    This is for extract P1 value from the matlab data with Josephson radiator's bias
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_std = []
        self.J_Bias = []
        self.J_Freq = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Projected_Occupation',
                             # data_type1='Weighted_Occupation',
                             data_type2='JB_Bias'):
        f0 = file_path[0]
        data0 = noiselib.loadmat(f0)
        occ_2D = data0[data_type1]
        J_Bias = data0[data_type2]
        J_Freq = J_Bias * 0.48
        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            if np.mean(occ_1D) >= 0.05:
                print('f=', f[-7: -4])
                # print('occ_1D=', occ_1D)
                print('occ_1D_avg=', np.mean(occ_1D))
            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters"""

        self.occ_1D_avg = np.average(occ_2D, axis=0)
        self.occ_1D_std = np.std(occ_2D, axis=0)
        self.J_Bias = np.array(J_Bias)
        self.J_Freq = np.array(J_Freq)


class P1_Avg_vs_Any(object):
    """
    This is for extract P1 value as function of any dependent
    """

    def __init__(self):
        self.p1_1D_avg = []
        self.p1_1D_std = []
        self.var = []
        self.var_name = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Weighted_Occupation',
                             data_type2='J1_Duration'):
        f0 = file_path[0]
        data0 = noiselib.loadmat(f0)
        p1_2D = data0[data_type1]
        var = data0[data_type2]
        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            p1_1D = np.array(data[data_type1])
            p1_2D = np.vstack((p1_2D, p1_1D))

        """Update parameters"""

        self.p1_1D_avg = np.average(p1_2D, axis=0)
        self.p1_1D_std = np.std(p1_2D, axis=0)
        self.var = np.array(var)
        self.var_name = data_type2

class P1_JSweep_Q2(object):
    """
    This is for extract P1 value from the matlab data with Josephson radiator's bias
    Q2 taken at different times
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_std = []
        self.J2_Bias = []
        self.J2_Freq = []

    def add_data_from_matlab(self, file_path1, file_path2,
                             data_type1='Projected_Occupation',
                             data_type2='J2_Bias'):

        f0 = file_path1[0]
        data0 = noiselib.loadmat(f0)
        occ_2D = data0[data_type1][18:]
        J2_Bias = data0[data_type2][18:]
        J2_Freq = J2_Bias * 0.48
        for f in file_path1[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1][18:])
            occ_2D = np.vstack((occ_2D, occ_1D))

        for f in file_path2:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters"""

        self.occ_1D_avg = np.average(occ_2D, axis=0)
        self.occ_1D_std = np.std(occ_2D, axis=0)
        self.J2_Bias = np.array(J2_Bias)
        self.J2_Freq = np.array(J2_Freq)


class GammaUp(object):
    """
    After with T1 and P1 values, get up transition rate
    """

    def __init__(self):
        self.temp = []
        self.QB_id = ''
        self.P1_2D = []  # [[P1(temp1), P1(temp1)], [P1(temp2), P1(temp2)], ...,]
        self.Gamma_Down_2D = []  # [[g(temp1), g(temp1)], [g(temp2), g(temp2)], ...,]
        self.Gamma_Up_2D = []
        self.Gamma_Up_1D = []
        self.Gamma_Up_1D_Error = []

    def add_data(self, temp, P1, Gamma_Down, QB_id):
        Gamma_Up = P1 * Gamma_Down / (1 - P1)
        self.temp.append(temp * 1.0 / 1000)
        self.P1_2D.append(P1)
        self.Gamma_Down_2D.append(Gamma_Down)
        self.Gamma_Up_2D.append(Gamma_Up)
        self.Gamma_Up_1D.append(np.mean(Gamma_Up))
        self.Gamma_Up_1D_Error.append(sem(Gamma_Up))
        self.QB_id = QB_id

    def _sort_temp(self):
        """
        sort data by temp if it is not imported
        :return:
        """
        temp = self.temp
        Gamma_Up_1D = self.Gamma_Up_1D
        Gamma_Up_1D_error = self.Gamma_Up_1D_Error

        zipped_lists_temp = zip(temp, temp)
        # print('zipped_lists_temp=', zipped_lists_temp)
        sorted_zipped_lists_temp = sorted(zipped_lists_temp)
        # print('sorted_zipped_lists_temp=', sorted_zipped_lists_temp)
        sorted_temp = np.array(
            [element for _, element in sorted_zipped_lists_temp])
        # print('sorted_temp=', sorted_temp)

        zipped_lists_Gamma_Up_1D = zip(temp, Gamma_Up_1D)
        sorted_zipped_lists_Gamma_Up_1D = sorted(zipped_lists_Gamma_Up_1D)
        # print('sorted_zipped_lists_Gamma_Up_1D=', sorted_zipped_lists_Gamma_Up_1D)
        sorted_Gamma_Up_1D = np.array(
            [element for _, element in sorted_zipped_lists_Gamma_Up_1D])
        # print('sorted_Gamma_Up_1D=', sorted_Gamma_Up_1D)

        zipped_lists_Gamma_Up_1D_error = zip(temp, Gamma_Up_1D_error)
        sorted_zipped_lists_Gamma_Up_1D_error = sorted(
            zipped_lists_Gamma_Up_1D_error)
        # print('sorted_zipped_lists_Gamma_Up_1D_error=', sorted_zipped_lists_Gamma_Up_1D_error)
        sorted_Gamma_Up_1D_error = np.array(
            [element for _, element in sorted_zipped_lists_Gamma_Up_1D_error])
        # print('sorted_Gamma_Up_1D_error=', sorted_Gamma_Up_1D_error)

        self.temp = sorted_temp
        self.GammaUp = sorted_Gamma_Up_1D
        self.Gamma_Up_1D_error = Gamma_Up_1D_error

    def plot(self):
        self._sort_temp()
        Temp = self.temp
        # print('type(Temp)=', type(Temp))
        QB_id = self.QB_id
        # print ('Temp=', Temp)
        Gamma_Up_1D = np.array(self.Gamma_Up_1D)
        Gamma_Up_1D_error = np.array(self.Gamma_Up_1D_Error)
        GammaErrorMinus = np.array(Gamma_Up_1D - Gamma_Up_1D_error)
        GammaErrorPlus = np.array(Gamma_Up_1D + Gamma_Up_1D_error)
        plt.scatter([1 / T for T in Temp],
                    [np.log(1 / G) for G in Gamma_Up_1D],
                    label='QB_Up')
        plt.fill_between([1 / T for T in Temp],
                         [np.log(1 / G) for G in GammaErrorMinus],
                         [np.log(1 / G) for G in GammaErrorPlus], alpha=0.2)
        # plt.plot([1 / T for T in Temp_array], [0.03 / T - 7.7 for T in Temp_array],
        #          label='Fit, f0=0.624GHz')

        plt.xlabel('1/Temp(Kelvin)')
        plt.ylabel('ln(1/Gamma)')
        plt.grid()
        plt.legend()
        plt.title(QB_id + ' Up Transition rate vs Temp')
        plt.show()

    def plot_iteration_at_temp(self):
        Temp = self.temp[0]
        QB_id = self.QB_id
        P1_2D = self.P1_2D[0]
        Gamma_Up_2D = self.Gamma_Up_2D[0]
        Gamma_Down_2D = self.Gamma_Down_2D[0]
        iteration = np.arange(0, len(P1_2D), 1)
        # print('P1_2D=', P1_2D)
        # print('Gamma_Up_2D=', Gamma_Up_2D)
        # print('Gamma_Down_2D=', Gamma_Down_2D)
        plt.plot(iteration, P1_2D / np.mean(P1_2D), label='P1')
        plt.plot(iteration, Gamma_Up_2D / np.mean(Gamma_Up_2D), label='Up')
        plt.plot(iteration, Gamma_Down_2D / np.mean(Gamma_Down_2D),
                 label='Down')
        plt.xlabel('Iteration')
        plt.ylabel('Normalized')
        plt.title(QB_id + ' at ' + str(Temp) + 'Kelvin')
        plt.grid()
        plt.legend()
        plt.ylim([0.5, 2])
        plt.show()


class TiedSphericalGaussianMixture(GaussianMixture):
    def _m_step(self, X, log_resp):
        super(TiedSphericalGaussianMixture, self)._m_step(X, log_resp)
        self.covariances_.fill(self.covariances_.mean())
        if np.any(np.less_equal(self.covariances_, 0.0)):
            raise ValueError(estimate_precision_error_message)
        self.precisions_cholesky_ = 1. / np.sqrt(self.covariances_)


class GMMFit(object):
    """
    For QB3 with 80 GHz mode, multiple higher states. Try to get a rate
    """

    def __init__(self):
        # self.statesCal = {
        #     "0Center": [0, 0],
        #     "0CenterStd": [0, 0],
        #     "1Center": [0, 0],
        #     "1CenterStd": [0, 0],
        #     "2Center": [0, 0],
        #     "2CenterStd": [0, 0],
        # }
        # self.time = []
        # self.IQ_I_Gate = []
        # self.IQ_X_Gate = []
        # self.IQ_pre = []
        # # 2D array [(time1)[I, Q], [I, Q], ...,
        # #           (time2)[I, Q], [I, Q], ...,
        # #           ...]
        # self.IQ = []    # same format as IQ_pre
        # self.IQ_filtered = [] # same format as IQ_pre

        self.I_Cal = {
            'iq_bucket': np.array([]),
            'gmm': None,
            'scale_factor': None,
            'means': None,
            'covs': None,
            'weights': None,
            'states': [0, 1, 2, 3, 4]
        }

        self.X_Cal = {
            'iq_bucket': np.array([]),
            'gmm': None,
            'scale_factor': None,
            'means': None,
            'covs': None,
            'weights': None,
            'states': [0, 1, 2, 3, 4]
        }

        self.J = {
            'iq_bucket': np.array([]),
            'gmm': None,
            'scale_factor': None,
            'means': None,
            'covs': None,
            'weights': None,
            'states': [0, 1, 2, 3, 4]
        }

        self.iq_bucket_I_Cal = None  # from Chris' fitting code
        self.iq_bucket_X_Cal = None  # from Chris' fitting code
        self.iq_bucket_J = None  # from Chris' fitting code

        # self.gmm = None
        # self.scale_factor = None
        # self.means = None
        # self.covs = None
        # self.weights = None

    def add_data_and_process_data(self, path):
        self._add_IQ_data_from_matlab(path=path)
        # self._fit_clusters(s='I')
        # self._fit_clusters(s='X')
        # self._fit_clusters(s='J')
        self._add_means_stds_weights_by_hand()

    def _add_IQ_data_from_matlab(self, path):
        data_I_Gate = noiselib.loadmat(path[0][0])
        data_X_Gate = noiselib.loadmat(path[1][0])
        data_J = noiselib.loadmat(path[2][0])

        IQ_I_Gate_I = data_I_Gate['Is']
        IQ_I_Gate_Q = data_I_Gate['Qs']
        IQ_X_Gate_I = data_X_Gate['Is']
        IQ_X_Gate_Q = data_X_Gate['Qs']
        IQ_J_I = data_J['Is']
        IQ_J_Q = data_J['Qs']
        iq_bucket_I_Cal = []
        iq_bucket_X_Cal = []
        iq_bucket_J = []

        for i in range(len(IQ_I_Gate_I)):
            iq_bucket_I_Cal.append([IQ_I_Gate_I[i], IQ_I_Gate_Q[i]])
            iq_bucket_X_Cal.append([IQ_X_Gate_I[i], IQ_X_Gate_Q[i]])
            iq_bucket_J.append([IQ_J_I[i], IQ_J_Q[i]])

        self.I_Cal['iq_bucket'] = np.array(iq_bucket_I_Cal)
        self.X_Cal['iq_bucket'] = np.array(iq_bucket_X_Cal)
        self.J['iq_bucket'] = np.array(iq_bucket_J)

    def _add_means_stds_weights_by_hand(self):
        means_I_Cal = np.array([
            [-4.98032745e-04, -5.29814968e-04],
            [-1.53501438e-05, -1.04117590e-03],
            [-1.18892628e-03, -7.60596566e-04],
            [-1.52127127e-05, -1.04145841e-03],
            [-1.29168464e-03, -1.36032532e-03]])
        covs_I_Cal = np.array([0.00018037, 0.00018037, 0.00018037, 0.00018037, 0.00018037])
        weights_I_Cal = np.array([0.20402284, 0.26886805, 0.0921508, 0.37787216, 0.05708614])

        means_X_Cal = np.array([
            [-8.12703158e-04, -5.26449187e-04],
            [-1.66068565e-05, -1.03175647e-03],
            [-1.27736318e-03, -1.40798388e-03],
            [-3.74014506e-04, -5.58047708e-04],
            [-1.25843911e-03, -8.65741341e-04]])
        covs_X_Cal = np.array([0.00016722, 0.00016722, 0.00016722, 0.00016722, 0.00016722])
        weights_X_Cal = np.array([0.107186, 0.2848422, 0.05143376, 0.46494838, 0.09158966])

        means_J = np.array([
            [-8.57133125e-04, -5.05998758e-04],
            [-2.14617765e-05, -1.03643611e-03],
            [-1.29915411e-03, -1.35906093e-03],
            [-4.29009115e-04, -5.40244316e-04],
            [-1.26567318e-03, -8.88872815e-04]])

        covs_J = np.array([0.00016538, 0.00016538, 0.00016538, 0.00016538, 0.00016538])
        weights_J = np.array([0.66665995, 0.05284824, 0.07662717, 0.07785456, 0.12601009])

        self.I_Cal['means'] = means_I_Cal
        self.I_Cal['covs'] = covs_I_Cal
        self.I_Cal['weights'] = weights_I_Cal

        self.X_Cal['means'] = means_X_Cal
        self.X_Cal['covs'] = covs_X_Cal
        self.X_Cal['weights'] = weights_X_Cal

        self.J['means'] = means_J
        self.J['covs'] = covs_J
        self.J['weights'] = weights_J

    def _fit_clusters(self, N=5, s='X'):
        """Fits N gaussian clouds to the N IQ blobs."""

        state = None
        if s == 'I':
            state = self.I_Cal
        elif s == 'X':
            state = self.X_Cal
        elif s == 'J':
            state = self.J

        combined_bucket = state['iq_bucket']
        # print('combined_bucket=', combined_bucket)  # LIU
        scale_factor = np.abs(combined_bucket.max())
        combined_bucket = combined_bucket / scale_factor
        minit = None  # np.array([np.mean(iq_bucket['Ground%s'%chan],axis=0),np.mean(iq_bucket['Excited%s'%chan],axis=0)])

        # init_means = np.array([
        #     [-3.74013889e-04, -5.58047925e-04],
        #     [-1.25843829e-03, -8.65740035e-04],
        #     [-1.66068104e-05, -1.03175658e-03],
        #     [-8.12700137e-04, -5.26448428e-04],
        #     [-1.27736337e-03, -1.40798313e-03]
        # ])
        #
        # init_covariances = np.array([0.00016722, 0.00016722, 0.00016722, 0.00016722, 0.00016722])

        gmm = TiedSphericalGaussianMixture(
            n_components=N,
            covariance_type='spherical',
            max_iter=300, tol=5e-7,
            means_init=minit, n_init=5)\
            .fit(combined_bucket)

        state['gmm'] = gmm
        state['scale_factor'] = scale_factor
        state['means'] = scale_factor * gmm.means_
        state['covs'] = scale_factor * np.sqrt(gmm.covariances_)
        state['weights'] = gmm.weights_

        print('means=', state['means'])
        print('covs=', state['covs'])
        print('weights', state['weights'])

        if s == 'I':
            self.I_Cal = state
        elif s == 'X':
            self.X_Cal = state
        elif s == 'J':
            self.J = state

    # def _get_gaus_means(self):
    #     means = self.scale_factor * self.gmm.means_
    #     return means
    #
    # #     ordered_states = self._id_centers_with_state_data(means
    # #     return {ordered_states[i]: means[i] for i in range(len(means))}
    # #
    # def _get_gaus_std(self):
    #     cov = self.scale_factor * np.sqrt(self.gmm.covariances_)
    #     return cov

    #     ordered_states = self._id_centers_with_state_data(means)
    #     return {ordered_states[i]: self.scale_factor * np.sqrt(
    #         self.gmm.covariances_[i])
    #             for i in range(len(means))}

    def _getCircle(self, IQCenter, std):
        """Returns points for state circle."""
        n = 50
        rxy = np.zeros((2, n))
        r0 = IQCenter
        r = std
        rxy[0] = r0[0] + r * np.cos(np.linspace(0, 2 * np.pi, n))
        rxy[1] = r0[1] + r * np.sin(np.linspace(0, 2 * np.pi, n))
        return rxy

    def plot_data(self, s='X', show=False):
        state = None
        if s == 'I':
            state = self.I_Cal
        elif s == 'X':
            state = self.X_Cal
        elif s == 'J':
            state = self.J

        IQCenters = state['means']
        stds = state['covs']
        I = state['iq_bucket'][:, 0]
        Q = state['iq_bucket'][:, 1]
        plt.figure(figsize=(7, 7))
        plt.scatter(I, Q, label='{}'.format(s), s=1)
        for i in range(len(IQCenters)):
            rxy = self._getCircle(IQCenters[i], stds[i])
            plt.scatter(rxy[0], rxy[1], marker='.')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.xlim(min(I) - 0.1 * (max(I) - min(I)),
                 max(I) + 0.1 * (max(I) - min(I)))
        plt.ylim(min(Q) - 0.1 * (max(Q) - min(Q)),
                 max(Q) + 0.1 * (max(Q) - min(Q)))
        plt.grid()
        plt.legend()
        plt.draw()
        if show:
            plt.show()


class GammaUp_New(object):
    """
    After with T1 and P1 values, get up transition rate. Avg P1/ T1 first then get up rate
    """

    def __init__(self):
        self.temp = []
        self.QB_id = ''
        self.P1_2D = []  # [[P1(temp1), P1(temp1)], [P1(temp2), P1(temp2)], ...,]
        self.Gamma_Down_2D = []  # [[g(temp1), g(temp1)], [g(temp2), g(temp2)], ...,]
        self.Gamma_Up_1D = []
        self.Gamma_Up_1D_Error = []

    def add_data(self, temp, P1, Gamma_Down, QB_id):
        print('here')
        P1_mean, P1_std = np.mean(P1), np.std(P1)
        Gamma_Down_mean, Gamma_Down_std = np.mean(Gamma_Down), np.std(
            Gamma_Down)
        Gamma_Up = P1_mean * np.mean(Gamma_Down_mean) / (1 - P1_mean)
        Gamma_Up_std = self.getGammaUpStd(P1_mean, P1_std, Gamma_Down_mean,
                                          Gamma_Down_std)
        self.temp.append(temp * 1.0 / 1000)
        self.Gamma_Up_1D.append(Gamma_Up)
        self.Gamma_Up_1D_Error.append(Gamma_Up_std)
        self.QB_id = QB_id

    def getGammaUpStd(self, P1_mean, P1_std, Gamma_Down_mean, Gamma_Down_std):
        x = (Gamma_Down_std / (1 / P1_mean - 1)) ** 2 + (
                Gamma_Down_mean * P1_std / ((P1_mean - 1) ** 2)) ** 2
        GammaUp_std = np.sqrt(x)
        # print('GammaUp_std=', GammaUp_std)
        return GammaUp_std

    def _sort_temp(self):
        """
        sort data by temp if it is not imported
        :return:
        """
        temp = self.temp
        Gamma_Up_1D = self.Gamma_Up_1D
        Gamma_Up_1D_error = self.Gamma_Up_1D_Error

        zipped_lists_temp = zip(temp, temp)
        # print('zipped_lists_temp=', zipped_lists_temp)
        sorted_zipped_lists_temp = sorted(zipped_lists_temp)
        # print('sorted_zipped_lists_temp=', sorted_zipped_lists_temp)
        sorted_temp = np.array(
            [element for _, element in sorted_zipped_lists_temp])
        # print('sorted_temp=', sorted_temp)

        zipped_lists_Gamma_Up_1D = zip(temp, Gamma_Up_1D)
        sorted_zipped_lists_Gamma_Up_1D = sorted(zipped_lists_Gamma_Up_1D)
        # print('sorted_zipped_lists_Gamma_Up_1D=', sorted_zipped_lists_Gamma_Up_1D)
        sorted_Gamma_Up_1D = np.array(
            [element for _, element in sorted_zipped_lists_Gamma_Up_1D])
        # print('sorted_Gamma_Up_1D=', sorted_Gamma_Up_1D)

        zipped_lists_Gamma_Up_1D_error = zip(temp, Gamma_Up_1D_error)
        sorted_zipped_lists_Gamma_Up_1D_error = sorted(
            zipped_lists_Gamma_Up_1D_error)
        # print('sorted_zipped_lists_Gamma_Up_1D_error=', sorted_zipped_lists_Gamma_Up_1D_error)
        sorted_Gamma_Up_1D_error = np.array(
            [element for _, element in sorted_zipped_lists_Gamma_Up_1D_error])
        # print('sorted_Gamma_Up_1D_error=', sorted_Gamma_Up_1D_error)

        self.temp = sorted_temp
        self.GammaUp = sorted_Gamma_Up_1D
        self.Gamma_Up_1D_error = Gamma_Up_1D_error

    def plot(self):
        self._sort_temp()
        Temp = self.temp
        QB_id = self.QB_id
        # print ('Temp=', Temp)
        Gamma_Up_1D = np.array(self.Gamma_Up_1D)
        Gamma_Up_1D_error = np.array(self.Gamma_Up_1D_Error)
        GammaErrorMinus = np.array(Gamma_Up_1D - Gamma_Up_1D_error)
        GammaErrorPlus = np.array(Gamma_Up_1D + Gamma_Up_1D_error)
        plt.scatter([1 / T for T in Temp],
                    [np.log(1 / G) for G in Gamma_Up_1D],
                    label='QB_Up')
        plt.fill_between([1 / T for T in Temp],
                         [np.log(1 / G) for G in GammaErrorMinus],
                         [np.log(1 / G) for G in GammaErrorPlus], alpha=0.2)
        plt.xlabel('1/Temp(Kelvin)')
        plt.ylabel('ln(1/Gamma)')
        plt.grid()
        plt.legend()
        plt.title(QB_id + ' Up Transition rate vs Temp')
        plt.show()

    def plot_iteration_at_temp(self):
        Temp = self.temp[0]
        QB_id = self.QB_id
        P1_2D = self.P1_2D[0]
        Gamma_Up_2D = self.Gamma_Up_2D[0]
        Gamma_Down_2D = self.Gamma_Down_2D[0]
        iteration = np.arange(0, len(P1_2D), 1)
        # print('P1_2D=', P1_2D)
        # print('Gamma_Up_2D=', Gamma_Up_2D)
        # print('Gamma_Down_2D=', Gamma_Down_2D)
        plt.plot(iteration, P1_2D / np.mean(P1_2D), label='P1')
        plt.plot(iteration, Gamma_Up_2D / np.mean(Gamma_Up_2D), label='Up')
        plt.plot(iteration, Gamma_Down_2D / np.mean(Gamma_Down_2D),
                 label='Down')
        plt.xlabel('Iteration')
        plt.ylabel('Normalized')
        plt.title(QB_id + ' at ' + str(Temp) + 'Kelvin')
        plt.grid()
        plt.legend()
        plt.ylim([0.5, 2])
        plt.show()


class calibration_debug(object):
    """
    Our P1 fluctuates with time. Is it from the measurement calibration
    Let's plot the IQ centers and std to see they fluctuate or not
    """

    def __init__(self):
        self.temp = None  # units mK
        self.gIavg = []
        self.gQavg = []
        self.gIstd = []
        self.gQstd = []
        self.eIavg = []
        self.eQavg = []
        self.eIstd = []
        self.eQstd = []

    def add_data_from_matlab(self, file_path, temp):
        gIavg = []
        gQavg = []
        gIstd = []
        gQstd = []
        eIavg = []
        eQavg = []
        eIstd = []
        eQstd = []
        for f in file_path:
            gIavg.append(noiselib.loadmat_ExptVars(f)['Ground_State_I'])
            gQavg.append(noiselib.loadmat_ExptVars(f)['Ground_State_Q'])
            gIstd.append(noiselib.loadmat_ExptVars(f)['Ground_State_Istd'])
            gQstd.append(noiselib.loadmat_ExptVars(f)['Ground_State_Qstd'])
            eIavg.append(noiselib.loadmat_ExptVars(f)['Excited_State_I'])
            eQavg.append(noiselib.loadmat_ExptVars(f)['Excited_State_Q'])
            eIstd.append(noiselib.loadmat_ExptVars(f)['Excited_State_Istd'])
            eQstd.append(noiselib.loadmat_ExptVars(f)['Excited_State_Qstd'])

        """Update parameters"""
        self.gIavg = np.array(gIavg)
        self.gQavg = np.array(gQavg)
        self.gIstd = np.array(gIstd)
        self.gQstd = np.array(gQstd)
        self.eIavg = np.array(eIavg)
        self.eQavg = np.array(eQavg)
        self.eIstd = np.array(eIstd)
        self.eQstd = np.array(eQstd)
        self.temp = temp


class BB_radiation(object):
    """
    To Understand more of the Blackbody radiation
    """

    def __init__(self):
        self.freq_range = np.linspace(1e8, 100e9, 500)
        self.temp_range = np.linspace(0.05, 0.5, 100)
        self.density_range = None
        self.temp = 3
        self.freq = 1e9

    def getDensity_freq(self):
        freq_range = self.freq_range
        print('k=', k)
        print('h=', h)
        temp = self.temp
        density_range = []
        for f in freq_range:
            density_range.append(
                (h * f) / (np.exp(h * f / (k * temp)) - 1))

        # density_range = np.array(
        #     (h * freq_range) / (np.exp(h * freq_range / (k * temp)) - 1))
        self.density_range = density_range

    def getDensity_temp(self):
        temp_range = self.temp_range
        freq = self.freq
        density_range = np.array(
            (h * freq) / (np.exp(h * freq / (k * temp_range)) - 1))
        self.density_range = density_range

    def plotDensity_freq(self):
        freq_range = self.freq_range
        density_range = self.density_range
        # plt.plot(freq_range, density_range)
        plt.loglog(freq_range, density_range)
        plt.xlabel('freq (Hz)')
        plt.ylabel('Scaled density')
        plt.show()

    def plotDensity_temp(self):
        temp_range = self.temp_range
        density_range = self.density_range
        plt.plot(temp_range, density_range / (1e9 * h),
                 label='Mode freq={} GHz'.format(self.freq / 1e9))
        plt.xlabel('Temp (Kelvin)')
        plt.ylabel('Scaled density')
        plt.legend()
        plt.show()


class Blackbox(object):
    """
    Understand the coupling between the 80 GHz mode to Qubit
    """

    def __init__(self):
        self.Antenna_mode = {
            "omega": None,
            "R": None,
            "C": None,
            "L": None,
            "Q": None,
            "T1": None
        }

        self.J = {
            "C": None,
            "L": None,
            "Ec": None
        }

        self.QB_mode = {
            "omega": None,
            "R": None,
            "C": None,
            "L": None,
            "Ec": None,
            "Q": None,
            "T1": None
        }
        self.Chi_QQ = None
        self.Chi_AA = None
        self.Chi_QA = None

        self.g = None

    def Input_update(self, AntennaInput, QBInput):
        self._Antenna_update(AntennaInput[0], AntennaInput[1], AntennaInput[2])
        self._QB_update(QBInput[0], QBInput[1], QBInput[2])
        self._getChi()
        self._getg()

    def _Antenna_update(self, f, C, R):
        """

        :param f: frequency
        :param C: effective capacitance; (1/2) ImY'(omega)
        :param R: effective resistance; 1/ReY(omega)
        :return: None; update the pole parameters
        """
        omega = 2 * pi * f
        L = 1 / (C * omega ** 2)
        Q = omega * C * R
        T1 = Q / omega
        self.Antenna_mode["omega"] = omega
        self.Antenna_mode["R"] = R
        self.Antenna_mode["C"] = C
        self.Antenna_mode["L"] = L
        self.Antenna_mode["Q"] = Q
        self.Antenna_mode["T1"] = T1

    def _QB_update(self, f, C, L):
        """

        :param f: frequency
        :param C: Total capcitance
        :param L: Josephson inductance
        :return: None; update the QB parameters
        """
        omega = 2 * pi * f
        C_J = 3.5 * 10 ** (-15)
        Ec_J = e ** 2 / (hbar * (2 * C_J))  # in units of omega
        Ec_QB = e ** 2 / (hbar * (2 * C))  # in units of omega
        self.QB_mode["omega"] = omega
        self.QB_mode["C"] = C
        self.QB_mode["L"] = L
        self.QB_mode["Ec"] = Ec_QB
        self.J["L"] = L
        self.J["C"] = C_J
        self.J["Ec"] = Ec_J

    def _getChi(self):
        Cj = self.J["C"]
        Lj = self.J["L"]
        Ec = self.J["Ec"]

        print('Ec=', Ec / (2 * pi * 1e9))

        C_Q = self.QB_mode["C"]
        L_Q = self.QB_mode["L"]
        C_A = self.Antenna_mode["C"]
        L_A = self.Antenna_mode["L"]

        Chi_QQ = -(L_Q / Lj) * (Cj / C_Q) * Ec
        Chi_AA = -(L_A / Lj) * (Cj / C_A) * Ec
        Chi_QA = -2 * np.sqrt(Chi_QQ * Chi_AA)

        self.Chi_QQ = Chi_QQ
        self.Chi_AA = Chi_AA
        self.Chi_QA = Chi_QA
        print("Chi_QQ=", Chi_QQ / (2 * pi * 1e9))
        print("Chi_AA=", Chi_AA / (2 * pi * 1e9))
        print("Chi_QA=", Chi_QA / (2 * pi * 1e9))

    def _getg(self):
        """
        We have Chi_QA = (g^2/Delta)*(alpha/(alpha+Delta))~= (g/Delta)^2 * alpha
        T_Q = (Delta/g)^2 * T_A Purcell limitation
        :return:
        """
        Chi_QA = self.Chi_QA
        Delta = self.Antenna_mode["omega"] - self.QB_mode["omega"]
        alpha = self.QB_mode["Ec"]
        T1_A = self.Antenna_mode["T1"]

        print("Delta=", Delta / (2 * pi * 1e9))
        print("alpha=", alpha / (2 * pi * 1e9))
        print("T1_A=", T1_A)

        factor = abs(Chi_QA / alpha)  # g^2/Delta^2
        g = np.sqrt(factor * Delta ** 2)
        T_Q = T1_A / factor

        self.QB_mode["T1"] = T_Q
        self.g = g


def getBBTensity(f, T):
    Intensity = (2 * h * f ** 3 / c ** 2) * (1 / (np.exp(h * f / (k * T)) - 1))
    return Intensity


def getGamma_pa(T, dfn=2e9, f0=120e9):
    """
    To calculate the theoretic photon assisted QP poisoning events based on
    the blackbody temperature and transmon's characteristic mode frequency
    :param T: temperature of the blackbody
    :param dfn: noise bandwidth
    :param f0: transmon's characteristic antenna mode frequency
    :return: Gamma_pa, the photon-assisted QP rate
    """
    Gamma_pa = dfn / (np.exp(h * f0 / (k * T)) - 1)
    return Gamma_pa