"""
PSD for Q1, Q2, Q4 at different temps
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
# from antennalib import getGamma_pa


def getGamma_pa(T, dfn=2e9, f0=120e9, Gammabase=1000):
    """
    To calculate the theoretic photon assisted QP poisoning events based on
    the blackbody temperature and transmon's characteristic mode frequency
    :param T: temperature of the blackbody
    :param dfn: noise bandwidth
    :param f0: transmon's characteristic antenna mode frequency
    :return: Gamma_pa, the photon-assisted QP rate
    """
    Gamma_pa = dfn / (np.exp(h * f0 / (k * T)) - 1)+Gammabase
    return Gamma_pa

Temp_QPT_Q1 = np.array([[0.076, 1429.9], [0.102, 1903.4], [0.199, 1897.7],
               [0.251, 2476.7], [0.284, 2922.8], [0.331, 5938.6],
               [0.372, 4121.9], [0.404, 6429], [0.451, 7073.7]])

pa_Q1 = getGamma_pa(Temp_QPT_Q1[:, 0], dfn=2e9, f0=92e9, Gammabase=1429)

Temp_QPT_Q2 = np.array([[0.076, 21.6], [0.102, 21.6], [0.199, 23.8],
                       [0.251, 25.5], [0.284, 31.1], [0.331, 98.9],
                       [0.372, 468.1], [0.404, 1305.2]])

pa_Q2 = getGamma_pa(Temp_QPT_Q2[:, 0], dfn=2e9, f0=120e9, Gammabase=21)

Temp_QPT_Q4 = np.array([[0.076, 326.6], [0.102, 335.0], [0.199, 375.9],
                       [0.251, 430.0], [0.284, 446.3], [0.331, 581.7],
                       [0.372, 1119.5], [0.404, 2773.2], [0.451, 8154.9]])
pa_Q4 = getGamma_pa(Temp_QPT_Q4[:, 0], dfn=2e9, f0=110e9, Gammabase=327)

Temp_QPT_Q4_old = np.array([[0.080, 1000.0/3.145], [0.185, 1000.0/2.847],
    [0.200, 1000.0/2.750], [0.210, 1000.0/2.709], [0.220, 1000.0/2.753],
    [0.229, 1000.0/2.683], [0.237, 1000.0/2.608], [0.245, 1000.0/2.649],
    [0.252, 1000.0/2.602], [0.260, 1000.0/2.441], [0.266, 1000.0/2.437],
    [0.274, 1000.0/2.469], [0.284, 1000.0/2.383], [0.294, 1000.0/2.228],
    [0.305, 1000.0/2.228], [0.318, 1000.0/2.165], [0.331, 1000.0/1.993],
    [0.353, 1000.0/1.526], [0.364, 1000.0/1.249], [0.375, 1000.0/1.049],
    [0.386, 1000.0/0.810], [0.396, 1000.0/0.654], [0.404, 1000.0/0.506],
    [0.412, 1000.0/0.422], [0.424, 1000.0/0.355]])

# plt.plot(1.0/Temp_QPT_Q1[:, 0], np.log(Temp_QPT_Q1[:, 1]), label='Q1')
# plt.plot(1.0/Temp_QPT_Q2[:, 0], np.log(Temp_QPT_Q2[:, 1]), label='Q2')
# plt.plot(1.0/Temp_QPT_Q4[:, 0], np.log(Temp_QPT_Q4[:, 1]), label='Q4')
plt.plot(Temp_QPT_Q1[:, 0], Temp_QPT_Q1[:, 1], label='Q1')
# plt.plot(Temp_QPT_Q1[:, 0], pa_Q1, '--', label='Q1_fit')
plt.plot(Temp_QPT_Q2[:, 0], Temp_QPT_Q2[:, 1], label='Q2')
# plt.plot(Temp_QPT_Q2[:, 0], pa_Q2, '--', label='Q2_fit')
plt.plot(Temp_QPT_Q4[:, 0], Temp_QPT_Q4[:, 1], label='Q4')
# plt.plot(Temp_QPT_Q4[:, 0], pa_Q4, '--', label='Q4_fit')
# plt.plot(Temp_QPT_Q4_old[:, 0], Temp_QPT_Q4_old[:, 1], 'k-', label='Q4_old')
plt.xlabel('Temp (Kelvin)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()
