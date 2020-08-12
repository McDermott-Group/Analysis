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

    def get_PSD(self, yale=True):
        parity_string_array = self.parity_string_array
        n = len(parity_string_array)
        for i in range(n):
            if yale:
                parity_string = (parity_string_array[i]-0.5)*2
                # parity_string = [(ps - 0.5) * 2 for ps in generate_hidden_signal(p_QP=[0.01, 0.01])]
            else:
                parity_string = parity_string_array[i]

            # f, Spp = welch(parity_string, fs) # welch seems to give weird result
            freq, spetral_power_density = periodogram(parity_string, self.fs)
            if i == 0:
                self.Spp_avg = spetral_power_density
                self.f_data = freq
            else:
                for j in range(len(self.Spp_avg)):
                    self.Spp_avg[j] += spetral_power_density[j]
        for i in range(len(self.Spp_avg)):
            self.Spp_avg[i] = self.Spp_avg[i]/n

    def plot_PSD(self, fit = False):
        self.get_PSD()
        fig = plt.figure()
        if fit:
            self.fit_PSD()
            plt.loglog(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0], self.params[1]), label='Fitted')
        plt.loglog(self.f_data, self.Spp_avg, label='Test')
        axes = plt.gca()
        axes.set_ylim([5e-6, 1e-2])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [1/Hz]')
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
# good data 07-27-20 On 600-1099, off 100-599

date = '08-12-20'
QP_Off_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/LIU/Q4_withQ5Poison/{}/QP_Tunneling_PSD_Poison_Off/MATLABData/')
QP_Off_files = np.arange(0, 199, 1)
QP_Off_filenames = [QP_Off_path.format(date) + 'QP_Tunneling_PSD_Poison_Off_{:03d}.mat'.format(i) for i in QP_Off_files]
QP_Off_psd = QP_PSD()
QP_Off_psd.add_data_from_matlab(QP_Off_filenames)
QP_Off_psd.plot_PSD(fit=False)


QP_On_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/LIU/Q4_withQ5Poison/{}/QP_Tunneling_PSD_Poison_On/MATLABData/')
QP_On_files = np.arange(0,199, 1) # 299 is good
QP_On_filenames = [QP_On_path.format(date) + 'QP_Tunneling_PSD_Poison_On_{:03d}.mat'.format(i) for i in QP_On_files]
QP_On_psd = QP_PSD()
QP_On_psd.add_data_from_matlab(QP_On_filenames)
QP_On_psd.plot_PSD(fit=False)

fig = plt.figure()
plt.loglog(QP_Off_psd.f_data, QP_Off_psd.Spp_avg, label='Poison Off')
plt.loglog(QP_On_psd.f_data, QP_On_psd.Spp_avg, label='Poison On')
axes = plt.gca()
axes.set_ylim([5e-5, 1e-2])
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [1/Hz]')
plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
plt.show()
# print QP_Off_psd.params