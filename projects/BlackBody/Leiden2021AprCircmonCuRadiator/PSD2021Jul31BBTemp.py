"""
PSD for Q1, Q2, Q4 at different temps
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.optimize import curve_fit
import copy
# from antennalib import getGamma_pa


def getGamma_pa(T, dfn, f0, Gammabase):
    """
    To calculate the theoretic photon assisted QP poisoning events based on
    the blackbody temperature and transmon's characteristic mode frequency
    :param T: temperature of the blackbody
    :param dfn: noise bandwidth
    :param f0: transmon's characteristic antenna mode frequency
    :return: Gamma_pa, the photon-assisted QP rate
    """
    # Gamma_pa = dfn / (np.exp(h * f0 / (k * T)) - 1)+Gammabase
    Gamma_pa = dfn * np.exp(-h * f0 / (k * T))+Gammabase
    return Gamma_pa

def fitToModel(Temp_QPT):
    """"""
    temp = Temp_QPT[:, 0]
    psd = Temp_QPT[:, 1]
    print('temp=', type(temp))
    # print('psd=', psd)
    # bounds = [(1e9, 100e9, 10), (3e9, 150e9, 1000)]
    params, covariance = curve_fit(getGamma_pa, temp, psd)
    return params

Temp_QPT_Q1 = np.array([
    [0.077, 0.55671], [0.199, 0.48771], [0.252, 0.39973], [0.275, 0.37894],
    [0.300, 0.32547], [0.325, 0.30938], [0.352, 0.27049], [0.375, 0.21793],
    [0.396, 0.17224], [0.428, 0.13323]])

Temp_QPT_Q2 = np.array([
    [0.077, 45.66544], [0.199, 41.41174], [0.252, 38.48387], [0.275, 35.45249],
    [0.300, 23.94614], [0.325, 11.36017], [0.352, 3.68751], [0.375, 1.27172],
    [0.396, 0.49732], [0.428, 0.18732]])

Temp_QPT_Q4 = np.array([
    [0.077, 2.80274], [0.199, 2.60836], [0.252, 2.28515], [0.275, 2.22227],
    [0.300, 2.07946], [0.325, 1.80078], [0.352, 1.24584], [0.375, 0.70726],
    [0.396, 0.38023], [0.428, 0.15798]])

Data_list = [Temp_QPT_Q1, Temp_QPT_Q2, Temp_QPT_Q4]
for d in Data_list:
    d[:, 1] = 1000/d[:, 1]

# params_Q2 = fitToModel(Temp_QPT_Q2)
# print('params_Q2=', params_Q2)
# pa_Q2 = getGamma_pa(Temp_QPT_Q2[:, 0], *params_Q2)

# pa_Q1 = getGamma_pa(Temp_QPT_Q1[:, 0], dfn=10000e9, f0=200e9, Gammabase=Temp_QPT_Q1[0][1])
pa_Q2 = getGamma_pa(Temp_QPT_Q2[:, 0], dfn=1e9, f0=120e9, Gammabase=Temp_QPT_Q2[0][1])
pa_Q4 = getGamma_pa(Temp_QPT_Q4[:, 0], dfn=10e9, f0=120e9, Gammabase=Temp_QPT_Q4[0][1])



# plt.plot(1.0/Temp_QPT_Q1[:, 0], np.log(Temp_QPT_Q1[:, 1]), label='Q1')
# plt.plot(1.0/Temp_QPT_Q2[:, 0], np.log(Temp_QPT_Q2[:, 1]), label='Q2')
# plt.plot(1.0/Temp_QPT_Q4[:, 0], np.log(Temp_QPT_Q4[:, 1]), label='Q4')
# plt.xlabel('1/Temp (Kelvin)')
# plt.ylabel('Ln(QPT (Hz))')
# plt.grid()
# plt.legend()
# plt.draw()
# plt.show()
plt.plot(Temp_QPT_Q1[:, 0], Temp_QPT_Q1[:, 1], label='Q1')
# plt.plot(Temp_QPT_Q1[:, 0], pa_Q1, '--', label='Q1_fit')
plt.plot(Temp_QPT_Q2[:, 0], Temp_QPT_Q2[:, 1], label='Q2')
plt.plot(Temp_QPT_Q2[:, 0], pa_Q2, '--', label='Q2_fit')
plt.plot(Temp_QPT_Q4[:, 0], Temp_QPT_Q4[:, 1], label='Q4')
plt.plot(Temp_QPT_Q4[:, 0], pa_Q4, '--', label='Q4_fit')
plt.xlabel('Temp (Kelvin)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.draw()
plt.show()

### the added part from Q3
# Extract the added channel
# Q1_added = copy.deepcopy(Temp_QPT_Q1)
# Q2_added = copy.deepcopy(Temp_QPT_Q2)
# Q4_added = copy.deepcopy(Temp_QPT_Q4)
# #
# Data_Added_list = [Q1_added, Q2_added, Q4_added]
# for d in Data_Added_list:
#     d[:, 1] = d[:, 1] - d[0][1]
# #
# # # print('Q2_20us_added=', Q2_20us_added)
# #
# # plt.figure()
# # plt.plot(Q1_added[1:, 0], Q1_added[1:, 1], label='Q1_added')
# # plt.plot(Q2_added[1:, 0], Q2_added[1:, 1], label='Q2_added')
# # plt.plot(Q4_added[1:, 0], Q4_added[1:, 1], label='Q4_added')
# # plt.xlabel('Temp (Kelvin)')
# # plt.ylabel('QPT (Hz)')
# # plt.yscale('log')
# # plt.grid()
# # plt.legend()
# # plt.show()
#
# plt.plot(1.0/Q1_added[1:, 0], np.log(Q1_added[1:, 0]), label='Q1')
# plt.xlabel('1/Temp (Kelvin)')
# plt.ylabel('Ln(QPT (Hz))')
# plt.grid()
# plt.legend()
# # plt.draw()
# plt.show()

