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
    def __init__(self, fs=1./(50e-6)):
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
                # parity_string = (parity_string_array[i]*0.5+0.5) # Chris's way
                # parity_string = generate_sin() # the sin function with
                # parity_string = [(ps - 0.5) * 2 for ps in generate_hidden_signal(p_QP=[0.05, 0.05])] # the simulated data works fine with the peridogram fft
                # parity_string = [ps for ps in generate_hidden_signal(p_QP=[0.05, 0.05])] # the simulated data works fine with the peridogram fft
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

    def plot_PSD(self, fit=False, plot=False):
        self.get_PSD()
        if plot:
            fig = plt.figure()
            if fit:
                self.fit_PSD()
                plt.loglog(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0], self.params[1]), label='Fitted')
                # plt.loglog(self.f_data, 2*self.fit_PSD_target_function(self.f_data, self.params[0]+120, self.params[1]*1.05), label='Fitted')
                # plt.plot(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0], self.params[1]), label='Fitted')
            plt.loglog(self.f_data, self.Spp_avg, label='Test fit [{:.2f}Hz], fidelity={:.2f}'.format(self.params[0], self.params[1]))
            # plt.plot(self.f_data, self.Spp_avg, label='Test fit [{:.2f}Hz], fidelity={:.2f}'.format(self.params[0], self.params[1]))
            axes = plt.gca()
            axes.set_ylim([1e-6, 1e-2])
            # axes.set_ylim([1e-6, 1e-1])
            plt.xlabel('frequency [Hz]')
            plt.ylabel('PSD [1/Hz]')
            plt.pause(0.1)
            plt.grid()
            plt.legend()
            plt.show()

    def fit_PSD_target_function(self, f, f_QP, f_RO):
        # return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) / self.fs
        return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2)

    def fit_PSD(self):
        # f_data = self.f_data
        # Spp_avg = self.Spp_avg
        # params = self.params
        self.f_data = self.f_data[1:]
        self.Spp_avg = self.Spp_avg[1:]

        initial_guess = [0.3, 0.5, 0.7]
        covariance = float('inf')
        for ig in initial_guess:
            # print f_data[:10]
            # sigma = self.f_data**1.8
            # sigma = self.Spp_avg**(-1)
            sigma = None
            params_curr, params_covariance_curr = optimize.curve_fit(
                self.fit_PSD_target_function, self.f_data, self.Spp_avg,
                bounds=[(10, 0), (10000, 1.2)], p0=[1000, ig], method='trf',
                sigma=sigma)
            print(params_covariance_curr)
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
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
# QP_files_Dirty = np.arange(10, 20, 1)
# DirtyAdd = np.arange(133, 136+1, 1)
# QP_files_Dirty = np.append(QP_files_Dirty, DirtyAdd)
# DirtyAdd2 = np.arange(167, 170+1, 1)
# DirtyAdd3 = np.arange(181, 182+1, 1)
# DirtyAdd4 = np.arange(236, 240+1, 1)
# DirtyAdd5 = np.arange(346, 349+1, 1)
# DirtyAdd6 = np.arange(408, 410+1, 1)
# DirtyAdd7 = np.arange(462, 465+1, 1)
# DirtyAdd8 = np.arange(534, 535+1, 1)
# DirtyAdd9 = np.arange(566, 570+1, 1)
#
# QP_files_Clean = np.arange(50, 60, 1)
# files_none = np.arange(1,n/11+1,1)*11-1
# QP_files = np.delete(QP_files, files_none)
# experiment_name = ('Interleave_PSD_Neg100')
# date = '10-28-20'
# QP_file_1 = np.arange(100, 110, 1)
# QP_file_1 = np.arange(23, 25, 1)
# file_1_add = np.arange(126, 130, 1)
# QP_file_1 = np.append(QP_file_1, file_1_add)
# file_1_add = np.arange(165, 170, 1)
# QP_file_1 = np.append(QP_file_1, file_1_add)
# file_1_add = np.arange(375, 380, 1)
# QP_file_1 = np.append(QP_file_1, file_1_add)
# file_1_add = np.arange(406, 410, 1)
# QP_file_1 = np.append(QP_file_1, file_1_add)
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_1]
# QP_psd_1 = QP_PSD()
# QP_psd_1.add_data_from_matlab(QP_filenames)
# QP_psd_1.plot_PSD(fit=True, plot=True)

# experiment_name = ('Interleave_PSD_Neg10')
# date = '10-29-20'
# QP_file_2 = np.arange(120, 130, 1)
# # QP_file_2 = np.arange(100, 120, 1)
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_2]
# QP_psd_2 = QP_PSD()
# QP_psd_2.add_data_from_matlab(QP_filenames)
# QP_psd_2.plot_PSD(fit=False)

# experiment_name = ('Interleave_PSD_Neg100')
# date = '10-28-20'
# QP_file_3 = np.arange(32, 35, 1)
# file_3_add = np.arange(52, 55, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# file_3_add = np.arange(121, 125, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# file_3_add = np.arange(147, 150, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# file_3_add = np.arange(192, 195, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# file_3_add = np.arange(201, 205, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# file_3_add = np.arange(461, 464, 1)
# QP_file_3 = np.append(QP_file_3, file_3_add)
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_3]
# QP_psd_3 = QP_PSD()
# QP_psd_3.add_data_from_matlab(QP_filenames)
# QP_psd_3.plot_PSD(fit=False)

# date = '10-27-20'
# QP_file_3 = np.arange(0, 100, 1)
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_3]
# QP_psd_3 = QP_PSD()
# QP_psd_3.add_data_from_matlab(QP_filenames)
# QP_psd_3.plot_PSD(fit=False)
#
# QP_file_4 = np.arange(500, 600, 1)
# QP_filenames = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_4]
# QP_psd_4 = QP_PSD()
# QP_psd_4.add_data_from_matlab(QP_filenames)
# QP_psd_4.plot_PSD(fit=False)

# QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}')
# date = '11-02-20'
# experiment_name = ('Interleave_PSD_Neg100')
# # QP_file_Clean = np.arange(442, 449, 1)     # 199 Good
# # QP_file_Clean = np.arange(70, 75, 1)     # 263 Good
# # QP_file_Clean = np.arange(340, 350, 1)     # 412 Hz, failed with Chris' PSD
# QP_file_Clean = np.arange(25, 30, 1)     # failed, failed with Chris' PSD
# QP_filenames_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file_Clean]
# QPT_Q6_NoPoison_Clean = QP_PSD()
# QPT_Q6_NoPoison_Clean.add_data_from_matlab(QP_filenames_Clean)
# QPT_Q6_NoPoison_Clean.plot_PSD(fit=True, plot=True)

# date = '11-04-20'
# experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(0, 10, 1)     # 221 Hz
# # QP_file = np.arange(0, 100, 1)     # 225 Hz
# # QP_file = np.arange(100, 150, 1)     # 234 Hz
# # QP_file = np.arange(210, 220, 1)     # 238 Hz
# QP_filenames_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_NoPoison_Clean = QP_PSD()
# QPT_Q6_NoPoison_Clean.add_data_from_matlab(QP_filenames_Clean)
# QPT_Q6_NoPoison_Clean.plot_PSD(fit=True, plot=True)

# experiment_name = ('Interleave_PSD_Neg10')
# QP_file = np.arange(0, 10, 1)     # 234 Hz
# # QP_file = np.arange(0, 50, 1)     # 229 Hz
# # QP_file = np.arange(0, 100, 1)     # 225 Hz
# # QP_file = np.arange(100, 150, 1)     # 234 Hz
# # QP_file = np.arange(210, 220, 1)     # 238 Hz
# QP_filenames_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
# QPT_Q6_NoPoison_Clean = QP_PSD()
# QPT_Q6_NoPoison_Clean.add_data_from_matlab(QP_filenames_Clean)
# QPT_Q6_NoPoison_Clean.plot_PSD(fit=True, plot=True)

date = '10-29-20'
QP_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')
experiment_name = ('Interleave_PSD_Neg100')
# QP_file = np.arange(10, 20, 1)     # 234 Hz
# QP_file = np.arange(440, 460, 1)     # 320 Hz
QP_file = np.arange(482, 490, 1)     # 320 Hz
QP_filenames_Clean = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in QP_file]
QPT_Q4_NoPoison_Clean = QP_PSD()
QPT_Q4_NoPoison_Clean.add_data_from_matlab(QP_filenames_Clean)
QPT_Q4_NoPoison_Clean.plot_PSD(fit=True, plot=True)



# def fit_PSD_target_function(f, f_QP, f_RO=1):
    # fs=1./100e-6
    # # return (1.0)*(4.0 * f_RO ** 2 * f_QP) / ((2.0* f_QP) ** 2 + (2.0 * np.pi * f) ** 2) + (1.0 - f_RO ** 2) / fs
    # return (4 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2)

    
# fig = plt.figure()
# plt.loglog(QP_psd_1.f_data, QP_psd_1.Spp_avg, label='No Poison')
# plt.loglog(QP_psd_2.f_data, QP_psd_2.Spp_avg, label='Poison ~ Al')
# # plt.loglog(QP_psd_3.f_data, QP_psd_3.Spp_avg, label='files with jump count rate~0.3')
# # plt.loglog(QP_psd_4.f_data, QP_psd_4.Spp_avg, label='4')
# print('sum(QP_psd_1.Spp_avg)=', sum(QP_psd_1.Spp_avg))
# print('sum(QP_psd_2.Spp_avg)=', sum(QP_psd_2.Spp_avg))
#
# axes = plt.gca()
# axes.set_ylim([1e-6, 1e-2])
# plt.xlabel('frequency [Hz]')
# plt.ylabel('PSD [1/Hz]')
# plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
# plt.show()
# # print QP_Off_psd.params