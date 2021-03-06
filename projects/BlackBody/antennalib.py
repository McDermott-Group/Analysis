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
            p1 = occ_2D[:, i]   # a list of Poplation at the same time [P1(ti), P2(ti), ...,Pn(ti)]
            avg, std = self._weighted_avg_and_std(p1, weights)
            occ_1D_avg.append(avg)
            error.append(std)
        # to get the error in the mean
        error = error/np.sqrt(len(weights))
        return np.array(occ_1D_avg), np.array(error)

    def _weighted_avg_and_std(self, values, weights):
        """
        Return the weighted average and standard deviation.

        values, weights -- Numpy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)
        return average, np.sqrt(variance)

    def fitToGammaUp(self):
        time = self.Pre_to_RO_Time
        P1 = self.occ_1D_avg
        if False:    # update length of the fit, use the first 25% of the data to fit
            length = int(len(time)/4)
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
            slope, intercept, r_value, p_value, std_err = linregress(time[:n], P1[:n])
            P1_fit = time * slope + intercept
            # print('P1_fit[:n]=', P1_fit[:n])
            for i in range(n):
                Chi_Square = Chi_Square + ((P1_fit[i]-P1[i])**2)/(P1_sig[i]**2)
            Chi_nu_Square = Chi_Square/(n-2)
            Chi_nu_Square_list.append(Chi_nu_Square)
            fit_number_list.append(n)
            GammaUp_list.append(slope)# units in Hz
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
        self.temp = np.append(self.temp, array_tempUpUpError[0])
        self.GammaUp = np.append(self.GammaUp, array_tempUpUpError[1])
        self.GammaUp_std_err = np.append(self.GammaUp_std_err, array_tempUpUpError[2])

    def plot(self):
        #
        Temp = self.temp * 10 ** (-3)
        Gamma = self.GammaUp
        GammaErrorMinus = self.GammaUp-self.GammaUp_std_err
        GammaErrorPlus = self.GammaUp+self.GammaUp_std_err
        # plt.scatter([1 / T for T in T_Pre],
        #             [np.log(1 / G) for G in Gamma_Pre_Q3], label='Q3_Up_Pre')
        # plt.plot([1/T for T in T_Pre_Q4], [np.log(1/G) for G in Gamma_Pre_Q4], label='Q4_Up_Pre')
        # plt.plot([1 / T for T in T_Pre], [0.2 / T - 8.7 for T in T_Pre],
        #          label='Fit, f0=4GHz')
        plt.plot([1 / T for T in Temp], [np.log(1 / G) for G in Gamma], label='Up_Pre')

        plt.fill_between([1 / T for T in Temp], [np.log(1 / G) for G in GammaErrorMinus],
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
        self.Gamma_fit_parameters = []   # save each file's fitting parameteres
        self.Gamma_fit_parameters = []   # save each file's fitting parameteres
        self.temp = None  # units mK
        self.fit_number_list = []
        self.T1_list = []

    def add_data_from_matlab(self, file_path, temp,
                             # data_type1='Weighted_Occupation',
                             data_type1='Projected_Occupation',
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
            Gamma_fit_parameters = np.hstack((Gamma_fit_parameters, fit_params))

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
        params, covariance = curve_fit(self.f_exp, time, P1, p0=[0.9, 30e4])
        # params, covariance = curve_fit(self.f_exp_off, time, P1, p0=[0.9, 30e4, 0.2])
        # print('params=', params)
        # print('covariance=', covariance)
        self.params = params
        return params

    def f_exp(self, t, amp, Gamma):
        return amp * np.exp(-t*Gamma)+0.09

    def f_exp_off(self, t, amp, Gamma, off):
        return amp * np.exp(-t*Gamma) + off

    def plot_histogram(self):
        Gamma = self.Gamma_fit_parameters
        T1 = [10**6/g for g in Gamma] # in us
        plt.hist(T1, bins = 10)
        plt.show()

    def plot(self):
        time = self.QB_Idle_Gate_Time * 10 ** 6
        occ = self.occ_1D_avg
        amp, Gamma = self.params[0], self.params[1]
        temp = str(self.temp) + 'mK'
        plt.figure(1)
        plt.plot(time, occ, 'k-', label=temp)
        # plt.plot(time, [amp*np.exp(-t*Gamma/(10**6)) for t in time], '--', label='fit')
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()


class P1(object):
    """
    This is for extract P1 value from the matlab data. Collaborate with T1
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.temp = None  # units mK

    def add_data_from_matlab(self, file_path, temp,
                             # data_type1='Weighted_Occupation'):
                             data_type1='Projected_Occupation'):
        occ_1D_avg = []
        for f in file_path:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_1D_avg.append(np.mean(occ_1D))

        """Update parameters"""
        # print('occ_1D_avg=', occ_1D_avg)

        self.occ_1D_avg = np.array(occ_1D_avg)
        self.temp = temp

class GammaUp(object):
    """
    After with T1 and P1 values, get up transition rate
    """
    def __init__(self):
        self.temp = []
        self.P1_2D = [] #[[P1(temp1), P1(temp1)], [P1(temp2), P1(temp2)], ...,]
        self.Gamma_Down_2D = [] #[[g(temp1), g(temp1)], [g(temp2), g(temp2)], ...,]
        self.Gamma_Up_2D = []
        self.Gamma_Up_1D = []
        self.Gamma_Up_1D_Error = []

    def add_data(self, temp, P1, Gamma_Down):
        Gamma_Up = P1*Gamma_Down/(1-P1)
        self.temp.append(temp)
        self.P1_2D.append(P1)
        self.Gamma_Down_2D.append(Gamma_Down)
        self.Gamma_Up_2D.append(Gamma_Up)
        self.Gamma_Up_1D.append(np.mean(Gamma_Up))
        self.Gamma_Up_1D_Error.append(sem(Gamma_Up))


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
