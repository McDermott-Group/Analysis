"""
A library for the antenna model
"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.stats import linregress
from scipy.stats import sem

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

    def add_data_from_matlab(self, file_path, temp,
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
        # plt.show(block=False)

    def plot(self):
        time = self.Pre_to_RO_Time * 10 ** 6
        occ = self.occ_1D_avg
        occ_fit = self.occ_1D_avg_fit
        error = self.occ_1D_avg_error_band
        GammaUp = self.GammaUp
        temp = str(self.temp) + 'mK'
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
        sorted_temp = np.array([element for _, element in sorted_zipped_lists_temp])
        print('sorted_temp=', sorted_temp)


        zipped_lists_GammaUp_std_err = zip(temp, GammaUp_std_err)
        sorted_zipped_lists_GammaUp_std_err = sorted(zipped_lists_GammaUp_std_err)
        sorted_temp_std_err = np.array([element for _, element in sorted_zipped_lists_GammaUp_std_err])

        zipped_lists_GammaUp = zip(temp, GammaUp)
        sorted_zipped_lists_GammaUp = sorted(zipped_lists_GammaUp)
        sorted_temp = np.array([element for _, element in sorted_zipped_lists_GammaUp])

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
                             data_type2='QB_Idle_Gate_Time'):
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

            fit_params = self.fitToExponential(occ_1D)[1]
            Gamma_fit_parameters = np.hstack(
                (Gamma_fit_parameters, fit_params))

            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters and units"""
        # occ_1D_avg, error = self._weighted_avg_and_std_2D(occ_2D, weights)
        # print('Gamma_fit_parameters=', Gamma_fit_parameters)

        self.occ_1D_avg = occ_1D
        self.temp = temp
        self.Gamma_fit_parameters = Gamma_fit_parameters

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

class P1_JSweep(object):
    """
    This is for extract P1 value from the matlab data with Josephson Radiator's bias
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.J2_Bias = []
        self.J2_Freq = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Projected_Occupation',
                             # data_type1='Phase',
                             # data_type1='Weighted_Occupation',
                             data_type2='J2_Bias'):

        occ_2D = []
        J2_Bias = []
        f0 = file_path[0]
        data0 = noiselib.loadmat(f0)
        occ_2D = data0[data_type1]
        # print('occ_2D=', occ_2D)
        J2_Bias = data0[data_type2]
        J2_Freq = J2_Bias*0.48
        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters"""
        # print('occ_1D_avg=', occ_1D_avg)

        self.occ_1D_avg = np.average(occ_2D, axis=0)
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
        sorted_zipped_lists_Gamma_Up_1D_error = sorted(zipped_lists_Gamma_Up_1D_error)
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
        sorted_zipped_lists_Gamma_Up_1D_error = sorted(zipped_lists_Gamma_Up_1D_error)
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


class BB_Radiation(object):
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
                 label='Mode freq={} GHz'.format(self.freq/1e9))
        plt.xlabel('Temp (Kelvin)')
        plt.ylabel('Scaled density')
        plt.legend()
        plt.show()

def getBBTensity(f, T):
    Intensity = (2*h*f**3/c**2)*(1/(np.exp(h*f/(k*T))-1))
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
