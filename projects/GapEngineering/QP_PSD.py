from noiselib import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram # there are also some other methods, chech them out
# from scipy.signal import welch
from numpy import random
from scipy import optimize

class QP_PSD(object):
    ###TO DO:
    #1. Fit
    def __init__(self, fs=1/100e-6):
        self.fs = fs
        self.parity_string_array = []
        self.Spp_avg = []
        self.f_data = None
        self.params = None

    def add_data_from_matlab(self, filenames, data_type='Single_Shot_Occupations'):
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

    def get_PSD(self):
        parity_string_array = self.parity_string_array
        for i in range(len(parity_string_array)):
            parity_string = parity_string_array[i]
            # parity_string = generate_hidden_signal(p_QP=[0.01, 0.01])
            # f, Spp = welch(parity_string, fs) # welch seems to give weird result
            freq, spetral_power_density = periodogram(parity_string, self.fs)
            if i == 0:
                self.Spp_avg = spetral_power_density
                self.f_data = freq
            else:
                for j in range(len(self.Spp_avg)):
                    self.Spp_avg[j] += spetral_power_density[j]

    def plot_PSD(self, fit = False):
        self.get_PSD()
        fig = plt.figure()
        if fit:
            self.fit_PSD()
            plt.loglog(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0]), label='Fitted')
        plt.loglog(self.f_data, self.Spp_avg, label='Test')
        axes = plt.gca()
        axes.set_ylim([1e-4, 1e-1])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [1/Hz]')
        plt.show()

    def fit_PSD_target_function(self, f, f_QP):
        f_RO = 0.999
        return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) / self.fs

    def fit_PSD(self):
        # f_data = self.f_data
        # Spp_avg = self.Spp_avg
        # params = self.params
        self.f_data = self.f_data[1:]
        self.Spp_avg = self.Spp_avg[1:]

        initial_guess = [100, 200]
        covariance = float('inf')
        for ig in initial_guess:
            # print f_data[:10]
            params_curr, params_covariance_curr = optimize.curve_fit(self.fit_PSD_target_function, self.f_data, self.Spp_avg, p0=[ig], method='trf')
            if params_covariance_curr[0] < covariance:
                self.params = params_curr
                params_covariance = params_covariance_curr
                covariance = params_covariance_curr[0]
        # print self.params


def generate_hidden_signal(length=8192, charge_burst_time=4000, p_QP=[0.01, 0.01]):
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


QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev01/Leiden_2020Feb/LIU/Q1/{}/QP_Tunneling_PSD_4/MATLABData/')
date = '03-16-20'
QP_files = np.arange(0, 8, 1)
filenames = [QP_path.format(date) + 'QP_Tunneling_PSD_4_{:03d}.mat'.format(i) for i in QP_files]

QP_psd = QP_PSD()
QP_psd.add_data_from_matlab(filenames)
QP_psd.plot_PSD()