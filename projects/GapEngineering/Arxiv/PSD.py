from scipy.signal import periodogram
from scipy import optimize
from numpy import random
import numpy as np
import matplotlib.pyplot as plt


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


def fit_PSD(f, f_QP):
    f_RO = 1
    return (4 * f_RO ** 2 * f_QP) / ((2 * f_QP) ** 2 + (2 * np.pi * f) ** 2) + (1 - f_RO ** 2) * 100e-6


Spp_avg = []
f_data = None
avg = 100
for i in range(avg):
    Hidden_Signal = generate_hidden_signal(p_QP=[0.01, 0.01])
    f, Spp = periodogram(Hidden_Signal, 10000)
    if i == 0:
        Spp_avg = Spp
        f_data = f
    else:
        for j in range(len(Spp_avg)):
            Spp_avg[j] += Spp[j]

# to get rid of the first 0
f_data = f_data[1:]
Spp_avg = Spp_avg[1:]

initial_guess = [100, 200]
covariance = float('inf')
params = None
for ig in initial_guess:
    params_curr, params_covariance_curr = optimize.curve_fit(fit_PSD, f_data, Spp_avg, p0=[ig], method='trf')
    if params_covariance_curr[0] < covariance:
        params = params_curr
        params_covariance = params_covariance_curr
        covariance = params_covariance_curr[0]


# print params

fig = plt.figure(figsize=(12, 4))
plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
plt.loglog(f_data, Spp_avg / avg, label='Simulated')
plt.loglog(f_data, fit_PSD(f_data, params[0]), label='Fitted')
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [1/Hz]')
plt.show()

