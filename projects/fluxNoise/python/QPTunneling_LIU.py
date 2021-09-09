"""
This is from Chris' QP Tunneling.py, This is just for test several things for Vincent Liu
"""

import os
import numpy as np
import noiselib
import matplotlib.pyplot as plt

from dataChest import *
from scipy.optimize import curve_fit
from scipy.signal import periodogram
#from Markov_Python2.analyze_QPTunneling_pomegranate import observed_to_recovered_signal_BW, generate_hidden_signal, hidden_to_observed_signal, generate_hidden_signal_highP1

reload(noiselib)


class QPTunneling_Wilen(object):

    def __init__(self, fs=1 / 50e-6, name='psd'):
        self.fs = fs
        vis = 1
        self.transfer = 1 / vis ** 2
        self.n_rows = 0
        self.name = name
        self.params = [0, 0]
        self.f_or_t = 'f'  # the PSD extraction is frequency or time
        self.T_parity = None
        self.one_over_f = False

    def add_datasets(self, file_path, data_str='Charge_Parity_Trace', simulate=False):
        if type(file_path) == str:
            file_path = [file_path]
        f_l = file_path[0]
        sample_rate = noiselib.loadmat_ExptVars(f_l)['Total_Time']
        sample_rate = sample_rate * 10**(-6)
        self.fs = 1.0/sample_rate
        for f in file_path:
            data = noiselib.loadmat(f)
            o = np.array(data[data_str])
            Serniak = True # This is used for fast measurement, M1 and M2's product
            # convertto10 = False   # This is for 1, -1 data to 0, 1 value conversion
            convertto10 = True
            if convertto10:
                o = o * 0.5 + 0.5
                for j in range(len(o)):
                    o[j] = map(int, o[j])
                # o = o * 0.5 + 0.5
            if Serniak:
                for j in range(len(o)):
                    innerProduct = [(2*o[j][k-1]-1)*(2*o[j][k]-1) for k in range(len(o[j]))]
                    # innerProduct = np.array(innerProduct)
                    innerProduct = map(int, innerProduct)
                    # print('type(innerProduct)=', type(innerProduct))
                    # o[j] = map(int, innerProduct)
                    o[j] = innerProduct
                # print('len(o[j])-1=', len(o[j])-1)
                # print('len(o[0])=', len(o[0]))
            if simulate:
                T_parity = 2 * 10 ** (-3)  # parity switching time
                p_QP = 1 - np.exp(-sample_rate / T_parity)  # converted to Poisson probability
                signal = generate_hidden_signal_highP1(p_QP=p_QP, P1=0.2)
                signal = hidden_to_observed_signal(signal, [0.9, 0.9])
                for i in range(len(o)):
                    o[i]=signal
            # print('o=', len(o[1]))
            self.add_data(o)

    def add_data(self, o):
        self.n_rows += 1  # don't divide by number of rows because these are already averaged in partition_and_avg_psd
        cpsd, f = noiselib.partition_and_avg_psd(o, self.fs)
        if not hasattr(self, 'cpsd'):
            self.cpsd, self.f = cpsd, f
        else:
            self.cpsd += cpsd

    def get_psd(self, window_averaging=True):
        cpsd = self.transfer * self.cpsd / self.n_rows
        if window_averaging:
            cpsd = noiselib.window_averaging(cpsd)
        return cpsd, self.f

    def get_fit_old(self, window_averaging=True):
        def y(x, a, gamma):
            return gamma * a / (np.pi ** 2 * x ** 2 + gamma ** 2)

        psd, f = self.get_psd(window_averaging)
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]
        popt, pcov = curve_fit(y, f, psd, bounds=(0, np.inf))
        # f_QP = self.params[0]
        # f_RO = self.params[1]
        # return fit_PSD_target_function(f, f_QP, f_RO), f
        self.params[0] = popt[-1]
        self.params[1] = 1
        return y(f, *popt), f

    def get_fit(self, window_averaging=True):
        def y(x, gamma, a):
            return gamma * a / (np.pi ** 2 * x ** 2 + gamma ** 2)

        psd, f = self.get_psd(window_averaging)
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]
        # print('len(f)=', len(f))
        #
        l = 5
        r = 500
        # print('f[l]', f[l])
        # print('f[r]', f[r])
        # f = f[l:]
        # psd = psd[l:]


        # params, covariance = curve_fit(y, f[l:r], psd[l: r], p0=[2000, 0.1],
        params, covariance = curve_fit(y, f, psd, p0=[5000, 0.1],
                                       bounds=[(1e1, 0), (1e6, 1)])# method='trf')

        print('params=', params)
        self.params = params
        f_QP = self.params[0]
        f_RO = self.params[1]
        self.f_or_t = 'f'
        self.T_parity = 1/f_QP
        # f = f[l:r]  # for partial fit
        return y(f, f_QP, f_RO), f

    # def get_fit_one_over_f(self, window_averaging=True):
    #     def y(x, gamma, a):
    #         return gamma * a / (np.pi ** 2 * x ** 2 + gamma ** 2)
    #
    #     psd, f = self.get_psd(window_averaging)
    #     psd = psd[~np.isnan(f)]
    #     f = f[~np.isnan(f)]
    #     f = f[~np.isnan(psd)]
    #     psd = psd[~np.isnan(psd)]
    #
    #     bounds = [(1e3, 0), (1e5, 1)]
    #
    #     # params, covariance = curve_fit(y, f, psd, p0=[0.1, 5000],
    #     #                                bounds=[(0, 1e3), (1, 1e5)], method='trf')
    #
    #     params, covariance = curve_fit(y, f, psd, bounds=bounds)
    #
    #     print('params=', params)
    #     self.params = params
    #     f_QP = self.params[0]
    #     f_RO = self.params[1]
    #     # amp = self.params[2]
    #     # alpha = self.params[3]
    #     self.f_or_t = 'f'
    #     self.T_parity = 1/f_QP
    #     # return y(f, 0, 0.0, 1e-4, 4e-1), f
    #     # return y(f, 6.55e3, 5.68e-2, amp, alpha), f
    #     return y(f, f_QP, f_RO), f

class QPTunneling_Liu(object):
    """
    This is done by periodogram method since this is actually conserving the integration (or the area)
    """

    def __init__(self, fs=1 / 50e-6, name='psd'):
        self.fs = fs
        self.parity_string_array = []
        self.psd_avg = []
        self.f_data = None
        self.params = None
        self.name = name
        self.f_or_t = 't'   # the PSD extraction is frequency or time
        self.T_parity = None

    def add_datasets(self, file_path, data_type='Charge_Parity_Trace', HMM=False, simulate=False):
        """

        :param file_path: a list of file names
        :param data_type:
        :return: array parity_string_arrays
        """
        f_l = file_path[0]
        sample_rate = noiselib.loadmat_ExptVars(f_l)['Total_Time']
        sample_rate = sample_rate * 10**(-6)
        # print('sample_rate=', sample_rate)
        # sample_rate = 50 * 10**(-6)
        self.fs = 1.0/sample_rate

        for f in file_path:
            data = noiselib.loadmat(f)
            ps_list = np.array(data[data_type])
            for ps in ps_list:
                if simulate:
                    T_parity = 2*10**(-3)   # parity switching time
                    p_QP = 1 - np.exp(-sample_rate/T_parity) # converted to Poisson probability
                    # signal = generate_hidden_signal(p_QP=[p_QP, p_QP])
                    signal = generate_hidden_signal_highP1(p_QP=p_QP)
                    # signal = generate_hidden_signal_two_freq()
                    signal = hidden_to_observed_signal(signal, [0.9, 0.9])
                    ps = [(ps - 0.5) * 2 for ps in signal]
                if HMM:
                    ps = self.apply_HMM(ps)
                self.parity_string_array.append(ps)

    def apply_HMM(self, parity_string):
        parity_string = [int(i * 0.5 + 0.5) for i in parity_string]
        parity_string = observed_to_recovered_signal_BW(parity_string)[0]
        parity_string = [int((i - 0.5) * 2) for i in parity_string]
        return parity_string

    def get_psd(self, window_averaging=False):
        parity_string_array = self.parity_string_array
        n = len(parity_string_array)

        for i in range(n):
            parity_string = parity_string_array[i]

            freq, psd = periodogram(parity_string, self.fs, return_onesided=True)
            if i == 0:
                self.psd_avg = psd
                self.f_data = freq
            else:
                for j in range(len(self.psd_avg)):
                    self.psd_avg[j] += psd[j]
        for i in range(len(self.psd_avg)):
            self.psd_avg[i] = self.psd_avg[i] / n

        single_to_double = 0.5  # FFT single to double sided amplitude
        self.psd_avg = single_to_double * self.psd_avg[1:]
        self.f_data = self.f_data[1:]

        if window_averaging:
            # print('self.psd_avg=', self.psd_avg)
            self.psd_avg = noiselib.window_averaging(self.psd_avg)
            # N = 20
            # self.psd_avg = np.convolve(self.psd_avg, np.ones((N,)) / N, mode='same')
        return self.psd_avg, self.f_data

    def get_fit(self):
        """
        Exactly the same as Serniak's thesis fit and parameters selection
        :return:
        """
        def fit_PSD_target_function(f, T_parity, F_map):
            return (4 * F_map ** 2 / T_parity) / ((2 / T_parity) ** 2 + (2 * np.pi * f) ** 2) + (1 - F_map ** 2) / self.fs

        psd, f = self.psd_avg, self.f_data
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]

        initial_guess = [0.3, 0.5, 0.7, 0.9]
        # initial_guess = [0.7]
        covariance = float('inf')
        for ig in initial_guess:
            sigma = f**1
            params_curr, params_covariance_curr = curve_fit(
                fit_PSD_target_function, f, psd,
                bounds=[(0.5*10**(-3), 0), (5*10**(-3), 1.0)], p0=[1.5*10**(-3), ig], method='trf',
                sigma=sigma)
            if params_covariance_curr[0][0] < covariance:
                self.params = params_curr
                params_covariance = params_covariance_curr
                covariance = params_covariance_curr[0][0]
        T_parity = self.params[0]
        F_map = self.params[1]
        self.T_parity = T_parity
        return fit_PSD_target_function(f, T_parity, F_map), f

    def get_fit_old(self):
        """
        Uses the frequency as the fit parameters
        :return:
        """
        def fit_PSD_target_function(f, f_QP, f_RO):
            tail = True
            # tail = False
            if tail:
                return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) / self.fs
            else:
                return (4 * f_RO ** 2 * f_QP) / \
                       ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2)

        psd, f = self.psd_avg, self.f_data
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]

        initial_guess = [0.3, 0.5, 0.7, 0.9]
        covariance = float('inf')
        for ig in initial_guess:
            sigma = f**1
            params_curr, params_covariance_curr = curve_fit(
                fit_PSD_target_function, f, psd,
                bounds=[(100, 0), (10000, 1.0)], p0=[100, ig], method='trf',
                sigma=sigma)
            if params_covariance_curr[0][0] < covariance:
                self.params = params_curr
                params_covariance = params_covariance_curr
                covariance = params_covariance_curr[0][0]
        f_QP = self.params[0]
        f_RO = self.params[1]
        return fit_PSD_target_function(f, f_QP, f_RO), f

#Dave's work on the fitting.
class QPTunneling_Harrison(object):
    """
    This is done by periodogram method since this is actually conserving the integration (or the area)
    """

    def __init__(self, fs=1 / 50e-6, name='psd'):
        self.fs = fs
        self.parity_string_array = []
        self.psd_avg = []
        self.f_data = None
        self.params = []
        self.name = name
        self.f_or_t = 't'   # the PSD extraction is frequency or time
        self.T_parity = []
        self.fidelity = []

    def add_datasets(self, file_path, data_type='Charge_Parity_Trace', HMM=False, simulate=False):
        """

        :param file_path: a list of file names
        :param data_type:
        :return: array parity_string_arrays
        """
        f_l = file_path[0]
        sample_rate = noiselib.loadmat_ExptVars(f_l)['Total_Time']
        sample_rate = sample_rate * 10**(-6)
        # print('sample_rate=', sample_rate)
        # sample_rate = 50 * 10**(-6)
        self.fs = 1.0/sample_rate

        for f in file_path:
            data = noiselib.loadmat(f)
            ps_list = np.array(data[data_type])
            for ps in ps_list:
                if simulate:
                    T_parity = 2*10**(-3)   # parity switching time
                    p_QP = 1 - np.exp(-sample_rate/T_parity) # converted to Poisson probability
                    # signal = generate_hidden_signal(p_QP=[p_QP, p_QP])
                    signal = generate_hidden_signal_highP1(p_QP=p_QP)
                    # signal = generate_hidden_signal_two_freq()
                    signal = hidden_to_observed_signal(signal, [0.9, 0.9])
                    ps = [(ps - 0.5) * 2 for ps in signal]
                if HMM:
                    ps = self.apply_HMM(ps)
                self.parity_string_array.append(ps)

    def apply_HMM(self, parity_string):
        parity_string = [int(i * 0.5 + 0.5) for i in parity_string]
        parity_string = observed_to_recovered_signal_BW(parity_string)[0]
        parity_string = [int((i - 0.5) * 2) for i in parity_string]
        return parity_string

    def get_psd(self, window_averaging=False, number=1, concatenate_records=1):
        parity_string_array = self.parity_string_array
        numRecords = len(parity_string_array)
        print(numRecords)
        if concatenate_records < 1:
            #this means to DIVIDE each record into 1/concatenate_records parts
            temp = []
            for record in parity_string_array:
                for splitRecord in np.split(record, int(1/concatenate_records)):
                    temp.append(splitRecord)
            parity_string_array = temp
        elif concatenate_records > 1:
            #this means to COMBINE concatenate_records
            temp = []
            for i in range(0, numRecords, concatenate_records):
                temp.append([])
                for j in range(0, concatenate_records):
                    temp[len(temp)-1] = np.concatenate([temp[len(temp)-1], parity_string_array[i+j]])
            parity_string_array=temp
        numRecords = len(parity_string_array)
        print(numRecords)

        #now we have a resized parity_string_array.  Now, we can computer the PSD of each string.
        #we will then divide the array into 'number" parts and  take each part and average.
        #We will return an array (length 'number') containing all of the PSDs

        #this computes the psd of each record and puts it into psds
        psds = []
        firstLoop = True
        for parity_string in parity_string_array:
            freq, psd = periodogram(parity_string, self.fs, return_onesided=True)
            psds.append(psd)
            if(firstLoop):
                firstLoop = False
                self.f_data = freq
                print(len(self.f_data))


        #this averages the proper number of psds and puts the averaged values into psd_ave
        psd_avg = []
        for i in np.arange(0, numRecords, int(1.0*numRecords/number)):
            psd_avg.append(np.zeros(len(psds[0])))
            for j in range(0, int(1.0*numRecords/number)):
                psd_avg[len(psd_avg) - 1] = np.add(psd_avg[len(psd_avg) - 1], psds[i + j])
            psd_avg[len(psd_avg) - 1] = psd_avg[len(psd_avg) - 1]/int(1.0*numRecords/number)

        self.psd_avg = psd_avg


        single_to_double = 0.5  # FFT single to double sided amplitude
        for i in range(0, len(psd_avg)):
            self.psd_avg[i] = single_to_double * self.psd_avg[i][1:]
            if window_averaging:
                # print('self.psd_avg=', self.psd_avg)
                self.psd_avg[i] = noiselib.window_averaging(self.psd_avg[i])
                # N = 20
                # self.psd_avg = np.convolve(self.psd_avg, np.ones((N,)) / N, mode='same')

        self.f_data = self.f_data[1:]
        return self.psd_avg, self.f_data

    def get_fit(self,excluded_points=1,ignore_fidelity=False):
        """
        Exactly the same as Serniak's thesis fit and parameters selection
        :return:
        """
        def fit_PSD_target_function(f, T_parity, F_map):
            return (4 * 1 * F_map ** 2 / T_parity) / ((2 / T_parity) ** 2 + (2 * np.pi * f) ** 2) + 1 * (1 - F_map ** 2) / self.fs
        self.T_parity=[]
        self.fidelity=[]
        toReturn = []

        f=[]
        for psd in self.psd_avg:
            f=self.f_data
            psd = psd[~np.isnan(f)]
            f = f[~np.isnan(f)]
            f = f[~np.isnan(psd)]
            psd = psd[~np.isnan(psd)]

            f=f[excluded_points:len(f)]
            psd=psd[excluded_points:len(psd)]
            initial_guess_knee = np.logspace(-5,-1,50)
            initial_guess_fidelity = np.linspace(0.2,1,8)
            # initial_guess = [0.7]
            rs = np.NINF #negative infinity is the lowest r2 value
            best_params=None
            for ig_knee in initial_guess_knee:
                for ig_fidelity in initial_guess_fidelity:
                    sigma = f**1
                    params_curr, params_covariance_curr = curve_fit(
                        fit_PSD_target_function, f, psd,
                        bounds=[(10**(-5), 0), (10**(-1), 1.0)], p0=[ig_knee, ig_fidelity],
                        sigma=sigma)
                    #just minimize r^2 to determine best fit. previously used covariance[0][0] which gave worse fits.
                    residuals = psd - fit_PSD_target_function(f, *params_curr)
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((psd - np.mean(psd)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot)
                    if rs < r_squared:
                    # if params_covariance_curr[0][0] < covariance:
                        rs = r_squared
                        best_params=params_curr
                        print(r_squared)
                        print(params_curr)
            #need to convert these to 2d
            self.params.append(best_params)
            self.T_parity.append(self.params[len(self.params)-1][0])
            self.fidelity.append(self.params[len(self.params)-1][1])
            toReturn.append(fit_PSD_target_function(f, self.T_parity[len(self.T_parity)-1], self.fidelity[len(self.fidelity)-1]))

        return toReturn, f

    def get_fit_old(self):
        """
        Uses the frequency as the fit parameters
        :return:
        """
        def fit_PSD_target_function(f, f_QP, f_RO):
            tail = True
            # tail = False
            if tail:
                return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) / self.fs
            else:
                return (4 * f_RO ** 2 * f_QP) / \
                       ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2)

        psd, f = self.psd_avg, self.f_data
        psd = psd[~np.isnan(f)]
        f = f[~np.isnan(f)]
        f = f[~np.isnan(psd)]
        psd = psd[~np.isnan(psd)]

        initial_guess = [0.3, 0.5, 0.7, 0.9]
        covariance = float('inf')
        for ig in initial_guess:
            sigma = f**1
            params_curr, params_covariance_curr = curve_fit(
                fit_PSD_target_function, f, psd,
                bounds=[(100, 0), (10000, 1.0)], p0=[100, ig], method='trf',
                sigma=sigma)
            if params_covariance_curr[0][0] < covariance:
                self.params = params_curr
                params_covariance = params_covariance_curr
                covariance = params_covariance_curr[0][0]
        f_QP = self.params[0]
        f_RO = self.params[1]
        return fit_PSD_target_function(f, f_QP, f_RO), f

class OneStateCleanDirty(object):
    """
    This is for clean dirty state file output
    """

    def __init__(self, fs=1 / 50e-6, name='psd'):
        self.fs = fs
        self.parity_string_array = []
        self.clean_files = []
        self.dirty_files = []
        self.clean_P1 = None
        self.dirty_P1 = None
        self.P1_threshold = 0.01
        self.P1_name = ""
        self.PSD_name = ""

    def update_experiment_name(self, P1_name, PSD_name):
        self.P1_name = P1_name
        self.PSD_name = PSD_name

    def add_p1_data_from_matlab(self, file_path,
                                data_type='Weighted_Occupation', P1_threshold=0.05):
                                # data_type='Single_Shot_Occupation', P1_threshold=0.05):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: clean and dirty file list
        """

        self.P1_threshold = P1_threshold
        P1_clean_list = []
        P1_dirty_list = []

        for f in file_path:
            data = noiselib.loadmat(f)
            p1_list = np.array(data[data_type])
            p1_avg = np.mean(p1_list)
            f_PSD = f.replace(self.P1_name, self.PSD_name)
            if p1_avg < self.P1_threshold:
                P1_clean_list.append(p1_avg)
                self.clean_files += [f_PSD]
            else:
                P1_dirty_list.append(p1_avg)
                self.dirty_files += [f_PSD]

        self.clean_P1 = np.mean(P1_clean_list)
        self.dirty_P1 = np.mean(P1_dirty_list)
        # print('clean list=', P1_clean_list)
        # print('dirty list=', P1_dirty_list)
        print('clean avg=', self.clean_P1)
        print('dirty avg=', self.dirty_P1)

def plotMultiFittedPSD(QPT_List, one_over_f=False, save=False, name=''):
    """
    plot multi fitted PSD for comparison
    :param QPT_List: a list of QPT object, e.g. [QPT_NoPoison, QPT_Neg10dBmPoison]
    :return: plot
    """
    plt.figure(figsize=(8, 6))
    linestyle = ['--', '+', 'o', 'v', 's', 'p', '*', 'h', 'x', 'D']
    for i, QPT in enumerate(QPT_List):
        # wa = False
        print('QPT=', QPT)
        wa = True
        psd, f = QPT.get_psd(window_averaging=wa)
        plt.loglog(f, psd, linestyle[i], label=r"{} PSD".format(QPT.name))
        # fit = False
        fit = True
        if fit:
            if one_over_f:
                print('oneoverf')
                psd_fit, f_fit = QPT.get_fit_one_over_f()
                F_map = QPT.params[1]
                if QPT.f_or_t == 't':
                    T_parity = QPT.params[0]
                else:
                    T_parity = (1/QPT.params[0])
                print('T_parity=', T_parity)
                print(psd_fit)
                plt.loglog(f_fit, psd_fit, '-',
                           label='{} fit [{:.5f} ms], fidelity={:.2f}'.format(
                               QPT.name, (T_parity)*10**3, F_map))
            else:
                psd_fit, f_fit = QPT.get_fit()
                F_map = QPT.params[1]
                if QPT.f_or_t == 't':
                    T_parity = QPT.params[0]
                else:
                    T_parity = (1/QPT.params[0])
                print('T_parity=', T_parity)
                plt.loglog(f_fit, psd_fit, '-',
                           label='{} fit [{:.5f} ms], fidelity={:.2f}'.format(
                               QPT.name, (T_parity)*10**3, F_map))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('S (1/Hz)')
    plt.grid()
    plt.grid()
    plt.grid()
    plt.legend()
    if save:
        plt.savefig(name+'.png')
    else:
        plt.show()
    # plt.pause(1)
    # plt.draw()


def plotMultiFittedPSD_Harrison(QPT_List, one_over_f=False, save=False, name='', concatenateRecords=1, excludedPoints=1, number=25):
    """
    plot multi fitted PSD for comparison
    :param QPT_List: a list of QPT object, e.g. [QPT_NoPoison, QPT_Neg10dBmPoison]
    :return: plot
    """
    plt.figure(figsize=(8, 6))
    linestyle = ['--', '+', 'o', 'v', 's', 'p', '*', 'h', 'x', 'D']
    for i, QPT in enumerate(QPT_List):
        # wa = False
        print('QPT=', QPT)
        wa = True
        psd, f = QPT.get_psd(number=number,window_averaging=wa,concatenateRecords=concatenateRecords)
        plt.loglog(f, psd, linestyle[i], label=r"{} PSD".format(QPT.name))
        # fit = False
        fit = True
        if fit:
            if one_over_f:
                print('oneoverf')
                psd_fit, f_fit = QPT.get_fit_one_over_f()
                F_map = QPT.params[1]
                if QPT.f_or_t == 't':
                    T_parity = QPT.params[0]
                else:
                    T_parity = (1/QPT.params[0])
                print('T_parity=', T_parity)
                print(psd_fit)
                plt.loglog(f_fit, psd_fit, '-',
                           label='{} fit [{:.5f} ms], fidelity={:.2f}'.format(
                               QPT.name, (T_parity)*10**3, F_map))
            else:
                psd_fit, f_fit = QPT.get_fit(excludedPoints=excludedPoints)
                F_map = QPT.params[1]
                if QPT.f_or_t == 't':
                    T_parity = QPT.params[0]
                else:
                    T_parity = (1/QPT.params[0])
                print('T_parity=', T_parity)
                plt.loglog(f_fit, psd_fit, '-',
                           label='{} fit [{:.5f} ms], fidelity={:.2f}'.format(
                               QPT.name, (T_parity)*10**3, F_map))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('S (1/Hz)')
    plt.grid()
    plt.grid()
    plt.grid()
    plt.legend()
    if save:
        plt.savefig(name+'.png')
    else:
        plt.show()
    # plt.pause(1)
    # plt.draw()