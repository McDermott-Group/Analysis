from noiselib import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram # there are also some other methods, chech them out
from scipy.signal import welch # welch gives me over filtered result
from numpy import random
from scipy import optimize
from QPTunneling import *

class QP_PSD(object):
    ###TO DO:
    #1. Fit
    # def __init__(self, fs=1./(100./1000000)):
    def __init__(self, fs=1./(100e-6)):
        self.fs = fs
        self.parity_string_array = []
        self.Spp_avg = []
        self.f_data = None
        self.params = None

    def add_data_from_matlab(self, filenames, data_type='Charge_Parity_Trace'):
        """

        :param filenames: a list of file names
        :param data_type:
        :return: array parity_string_arrays
        """
        for file in filenames:
            data = loadmat(file)
            ps_list = np.array(data[data_type])
            for ps in ps_list:
                self.parity_string_array.append(ps)

    def get_PSD(self, yale=True):
        parity_string_array = self.parity_string_array
        n = len(parity_string_array)
        for i in range(n):
            if yale:
                # parity_string = (parity_string_array[i]*1.0-0.5)*2
                parity_string = parity_string_array[i] # Chris's way
                # parity_string = generate_sin() # the sin function with
                # parity_string = [(ps - 0.5) * 2 for ps in generate_hidden_signal(p_QP=[0.01, 0.01])] # the simulated data works fine with the peridogram fft
                # parity_string = [ps for ps in generate_hidden_signal(p_QP=[0.01, 0.01])] # the simulated data works fine with the peridogram fft
            else:
                parity_string = parity_string_array[i]

            # freq, spetral_power_density = welch(parity_string, self.fs) # welch seems to give weird result
            freq, spetral_power_density = periodogram(parity_string, self.fs) # this passes the sin
            # o = np.array(list([parity_string])) # Chris's way
            # spetral_power_density, freq = noiselib.partition_and_avg_psd(o, self.fs)
            if i == 0:
                self.Spp_avg = spetral_power_density
                self.f_data = freq
            else:
                for j in range(len(self.Spp_avg)):
                    self.Spp_avg[j] += spetral_power_density[j]
        for i in range(len(self.Spp_avg)):
            self.Spp_avg[i] = self.Spp_avg[i]/n

    def plot_PSD(self, fit = False, plot = False):
        self.get_PSD()
        if plot:
            fig = plt.figure()
            if fit:
                self.fit_PSD()
                plt.loglog(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0], self.params[1]), label='Fitted')
            plt.loglog(self.f_data, self.Spp_avg, label='Test')
            axes = plt.gca()
            axes.set_ylim([1e-6, 1e-2])
            # axes.set_ylim([1e-6, 1e-1])
            plt.xlabel('frequency [Hz]')
            plt.ylabel('PSD [1/Hz]')
            plt.pause(0.1)
            plt.show()

    def fit_PSD_target_function(self, f, f_QP, f_RO):
        return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) / self.fs

    def fit_PSD(self):
        # f_data = self.f_data
        # Spp_avg = self.Spp_avg
        # params = self.params
        self.f_data = self.f_data[1:]
        self.Spp_avg = self.Spp_avg[1:]

        initial_guess = [0.99, 0.999, 0.9999]
        covariance = float('inf')
        for ig in initial_guess:
            # print f_data[:10]
            params_curr, params_covariance_curr = optimize.curve_fit(self.fit_PSD_target_function, self.f_data, self.Spp_avg, p0=[100, ig], method='trf')
            print params_covariance_curr
            if params_covariance_curr[0][0] < covariance:
                self.params = params_curr
                params_covariance = params_covariance_curr
                covariance = params_covariance_curr[0][0]

def generate_sin(n = np.arange(0, 10000,1.0), fs=10000.0, f_QP = 200):
    # 10000 data points
    n = n/fs
    # 100 Hz
    s = np.sin(2*np.pi*f_QP*n)
    return s
def generate_hidden_signal(length=10000, charge_burst_time=4000, p_QP=[0.01, 0.01]):
    """
    Generate hidden signal from given parameters
    :param length: the length of the hidden signal, eg. 8192
    :param charge_burst_time: the time where we have the charge burst, eg. 2456
    :param p_QP: p_QP = [p_QP_before, p_QP_after], a two element list consists QP tunneling rate before and after the burst
        p_QP = 1.0-np.exp(-t_meas/t_QP), e.g. t_meas = 100 us, t_QP = 10 ms, p_QP=0.01
    :return: the hidden signal [0,1,0,0,1,1,1,...]
    """

    p_QP_before = p_QP[0]
    p_QP_after = p_QP[1]
    hidden_signal = [0] * length

    hidden_signal[0] = random.choice([0, 1], p=[0.5, 0.5])  # Initial parity

    for i in range(1, charge_burst_time):
        hidden_signal[i] = random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                         p=[1 - p_QP_before, p_QP_before])
    for i in range(charge_burst_time, length):
        hidden_signal[i] = random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                         p=[1 - p_QP_after, p_QP_after])

    return hidden_signal
# good data 07-27-20 On 600-1099, off 100-599
# n = 10
# date = '08-23-20'
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
QP_files = np.arange(110, 140, 1)
# files_none = np.arange(1,n/11+1,1)*11-1
# QP_files = np.delete(QP_files, files_none)

date = '10-18-20'
experiment_name = ('Interleave_PSD')
QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_files]
QP_psd = QP_PSD()
QP_psd.add_data_from_matlab(QP_filenames)
QP_psd.plot_PSD(fit=False, plot=True)


# def fit_PSD_target_function(f, f_QP, f_RO=1):
    # fs=1./100e-6
    # # return (1.0)*(4.0 * f_RO ** 2 * f_QP) / ((2.0* f_QP) ** 2 + (2.0 * np.pi * f) ** 2) + (1.0 - f_RO ** 2) / fs
    # return (4 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2)

    
# fig = plt.figure()
# plt.loglog(QP_psd.f_data, QP_psd.Spp_avg, label='Neg2 dBm')
#
# axes = plt.gca()
# axes.set_ylim([1e-6, 1e-2])
# plt.xlabel('frequency [Hz]')
# plt.ylabel('PSD [1/Hz]')
# plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
# plt.show()
# # print QP_Off_psd.params