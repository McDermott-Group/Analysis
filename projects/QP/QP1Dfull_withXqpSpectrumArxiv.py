#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:00:25 2021

@author: robertmcdermott

Modified:
2022Jan
Chuanhong Liu

For the up transition rate calculation

"""

# from pylab import plot, xlabel, ylabel, show
import matplotlib.pyplot as plt
import random
import time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
import cmath

Delta = 1  # energy gap
L = 2  # length in microns, assume in 2 um, the n_QP is uniform distribute
tau0 = 400  # e-p coupling time, units of ns
Dn = 6e-2  # units of um^2/ns -- use 6 here for aluminum

nx = 1  # points in space grid, let's say it is 1 D now
ne_log = 10  # points in log energy grid
# nt = 10000*500  # time steps

time_tot = 0.00001*1e-3   # milli sec. 1 msec should be enough, from 1 to 10 msec,
                    # the change is 0.001% level, 100 msec is more than enough without trapping
s = 1e-2  # sets time step; units of tau0 -- try 1e-4, 5e-2 good for E=2Delta injection
dt = s * tau0
Gamma = 1e-3  # Dynes broadening

nt = int(time_tot*1e9/dt)  # time steps, convert things to nanoseconds

# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
emax = 10.0  # units of Delta
# de = (emax - 1.0) / ne  # energy increment in units of Delta = 1
# e = 1.0 + de * np.linspace(0, ne - 1, ne)
e_log = np.logspace(0, 4, num=ne_log, base=2)
e_inj = 2
idx_e_upper = (np.abs(e_log-e_inj)).argmin()
e = e_log[:idx_e_upper+5] # useful energy grid
ne = len(e) # points in energy grid
broadened = complex(0, -Gamma)
Dynes = e + broadened
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real
n = np.zeros((nx, ne))
# n = np.zeros(ne)
# print('dimension=', n.ndim)

G = np.zeros((ne, ne))
Gr0 = np.zeros((ne, ne))
Gp = np.zeros((ne, ne))

inflow = np.zeros(ne)
outflow = np.zeros(ne)
phonons_s = np.zeros(ne)

# generate the scattering and recombination matrices

### kTc = \Delta/1.76
for i in range(ne):
    for j in range(ne):
        G[i, j] = s * (1.76 ** 3.0) * ((e[i] - e[j]) ** 2.0) * np.heaviside(e[j] - e[i], 0.5) * rho[j] \
                  * (1.0 - Delta ** 2.0 / (e[i] * e[j]))
        # this is scattering from j to i
        Gr0[i, j] = s * (1.76 ** 3.0) * ((e[i] + e[j]) ** 2.0) * (1.0 + Delta ** 2.0 / (e[i] * e[j]))

# generate the scattering matrix for phonons
for i in range(ne):
    for j in range(ne):
        if i <= j:
            Gp[i, j] = G[j - i, j]

logG = np.log(G)
GT = G.transpose()
logGT = np.log(GT)
# print(GT)
column = np.ones(ne)
column.shape = (ne, 1)
outflow_helper = np.matmul(GT, column)
# print(outflow_helper)

def scatter(n, inj, dt):
    """
    It takes the recombination process into account
    :param n:
    :param inj:
    :param dt:
    :return:
    """
    print('scatter')
    print('n=', n)
    n.shape = (len(n), 1)
    print('After shape1')
    print('n=', n)
    inflow = np.matmul(G, n)  # this is scattering into the state
    # print('inflow=', inflow)
    outflow = outflow_helper * n  # elementwise multiplication
    # print('outflow=', outflow)
    recomb_helper = np.kron(n, np.ones(ne))
    Gr = 2 * Gr0 * recomb_helper    # 2 QPs recombine into 1 CP
    # print('Gr=', Gr)
    # print('n=', n)
    # print('here')
    recombined = np.matmul(Gr, n)
    # print('recombined=', recombined)
    n = n + inflow - outflow - recombined
    n = n + inj
    print('After dynamics')
    print('n=', n)
    print('len(n)=', len(n))
    n.shape = (1, len(n))
    print('After shape2')
    print('n=', n)
    return n

# injection
def inj(ne, e_inj, Gamma_p, tau0, dt):
    """

    :param ne:  energy bins
    :param e_inj:
    :param Gamma_p: the measured parity rate
    :return: inj, the injection array in units of Lenander paper
    """
    delta_inj = np.zeros(ne)
    delta_inj.shape = (ne, 1)
    idx = (np.abs(e-e_inj)).argmin()

    Vol_Al = 2.0  #um^3
    t = tau0*dt*1e-9    # in units of seconds
    N_qp = 2 * t * Gamma_p  # time * parity rate = QPs generated, 2 means a broken pair generates 2 QPs
    # print('N_qp=', N_qp)
    n_cp = 4.0*1e6    # Cooper pair density um^-3
    n_qp = N_qp/(Vol_Al)
    n_inj = (Delta/2.0)*(n_qp/n_cp)
    delta_inj[idx] = n_inj # number of QPs, not the reduced x_QP, this only linear scales the result

    return delta_inj

Gamma_p = 4.4e3
n_inj = inj(ne, e_inj, Gamma_p, tau0, dt)

### evolve in time
for i in range(nt):
    if i % 10000 == 0:
        print(i/10000)
        print('n[0, :]=', n[0, :])
    # if i % 1000 == 0:
    #     print(i/1000)
    # if i % 1000 == 0:
    #     plt.plot(e, n[0, :])
    #     plt.xlabel('$E/\Delta$')
    #     plt.ylabel('QP Density')
        # plt.yscale('log')
        # plt.show()
    # no spatial consideration here
    inj = n_inj
    n[0, :] = scatter(n[0, :], inj, dt)

# print('sum(n[0, :])=', sum(n[0, :]))
print('x_qp=', 2*sum(n[0, :]))
ni = n[0, :]

QP_Data = []
for i in range(len(ni)):
    d = [e[i], ni[i]]
    QP_Data.append(d)
np.savetxt('QPEnergySpectrum_250GHz.txt', QP_Data)

plt.plot(e, ni)
plt.xlabel('$E/\Delta$')
plt.ylabel('ni')
plt.yscale('log')
plt.show()



# e_phon = e - 1.0
# e_recomb = 2.0 + np.linspace(0, 2 * ne - 2, 2 * ne - 1)
#
# # look at spectrum of scattered phonons
# phonons_s = np.matmul(Gp, n[nx // 2, :])
# phonons_r = np.zeros(2 * ne - 1)
#
# # look at spectrum of recombination phonons
# for k in range(2 * ne - 1):
#     l = 0
#     while (l <= k // 2) and (k // 2 + l < ne - 1):  # need to fix this -- treat even and odd cases separately
#         # think this is not quite right -- indexing probably screwed up slightly at top of energy range
#         if k % 2 == 0:  # k even, sum from main diagonal
#             phonons_r[k] = phonons_r[k] + Gr0[k // 2 - l, k // 2 + l] * n[nx // 2, k // 2 - l] * n[nx // 2, k // 2 + l]
#         else:
#             phonons_r[k] = phonons_r[k] + Gr0[k // 2 - l, k // 2 + l + 1] * n[nx // 2, k // 2 - l] * n[
#                 nx // 2, k // 2 + l + 1]
#         l += 1
#
# e_phon = e - 1
# e_recomb = 2 + 2 * (emax - 1) * np.linspace(0, 1, 2 * ne - 1)
#
# # plot(e_phon, np.log(phonons_s), 'b', e_recomb, np.log(phonons_r), 'r')
# plt.plot(e_phon, phonons_s, 'b', label='scatter')
# plt.plot(e_recomb, phonons_r, 'r', label='recomb')
# plt.yscale('log')
# plt.legend()
# plt.show()
