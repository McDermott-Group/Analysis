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


def Chi_Squaire(n_curr, n_before):
    ChiSquare = 0.0
    # for i in range(len(n_curr) - 4):
    #     chisquare = (n_curr[i] - n_before[i]) ** 2.0 / n_curr[i]
    #     ChiSquare = ChiSquare + chisquare
    return ChiSquare

def myConverge(n_curr, n_before):
    myC = 1
    for i in range(len(n_curr) - 4):
        res = np.abs(n_curr[i] - n_before[i])/n_curr[i]
        if res > 1.0e-10:   # or 1e-9
            # print('res=', res)
            myC = 0
            break
    return myC


Delta = 1  # energy gap
L = 2  # length in microns, assume in 2 um, the n_QP is uniform distribute
tau0 = 400  # e-p coupling time, units of ns
Dn = 6e-2  # units of um^2/ns -- use 6 here for aluminum

nx = 1  # points in space grid, let's say it is 1 D now
ne_log = 100    # points in log energy grid,
                # 1000 enough for my purpose, 100 good for debug purpose
                # finer ne also means finer time step

time_tot = 0.01 * 1e-3  # milli sec. 1 msec should be enough, from 1 to 10 msec,
# the change is 0.001% level, 100 msec is more than enough without trapping
s = 1e-4  # sets time step; units of tau0 -- try 1e-4, 5e-2 good for E=2Delta injection
dt = s * tau0   # time step, units of nsec
Gamma = 1e-3  # Dynes broadening

nt = int(time_tot * 1e9 / dt)  # time steps, convert things to nanoseconds

# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
emax = 10.0  # units of Delta
# de = (emax - 1.0) / ne  # energy increment in units of Delta = 1
# e = 1.0 + de * np.linspace(0, ne - 1, ne)
e_log = np.logspace(0, 4, num=ne_log, base=2)
e_inj = 5.43/2
idx_e_upper = (np.abs(e_log - e_inj)).argmin()
e = e_log[:idx_e_upper + 5]  # useful energy grid
ne = len(e)  # points in energy grid
broadened = complex(0, -Gamma)
Dynes = e + broadened
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real
n = np.zeros(ne)


### initial data
dat = np.loadtxt('QPEnergySpectrum_250GHzPhoton.txt')
e = dat[:, 0]
n = dat[:, 1]# use earlier data as initial

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
        # this is scattering from j to i
        G[i, j] = s * (1.76 ** 3.0) * ((e[i] - e[j]) ** 2.0) * np.heaviside(e[j] - e[i], 0.5) * rho[j] \
                  * (1.0 - Delta ** 2.0 / (e[i] * e[j]))
        # recombination between i and j energy bins
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

def scatter(n, inj, dt, i):
    """
    It takes the recombination process into account
    :param n:
    :param inj:
    :param dt:
    :return:
    """
    n.shape = (len(n), 1)
    inflow = np.matmul(G, n)  # this is scattering into the state
    outflow = outflow_helper * n  # elementwise multiplication
    recomb_helper = np.kron(n, np.ones(ne))  # helper is the issue, get large numbers and become nans
    Gr = 2 * Gr0 * recomb_helper  # 2 QPs recombine into 1 CP
    recombined = np.matmul(Gr, n)

    trap_time = 0.01e-6  # 100 usec
    trap_rate = 1/(trap_time*1e9)   # convert to 1/ns
    trapped = dt*trap_rate*n
    # if i % 10000 == 0:
    #     print(i / 10000)
    #     print('sum(recombined)=', sum(recombined))
    #     print('sum(trapped)=', sum(trapped))

    # n = n + inflow - outflow - recombined-trapped
    n = n + inflow - outflow - recombined
    n = n + inj

    ### Debug
    # nn = 1
    # if i % nn == 0 and i <2:
    #     print(i/nn)
    #     # print('n=', n)
    #     print('sum(n)=', sum(n))
    #     # if sum(inflow) < 0:
    #     ### issue, some of the results are negative, which is not physical
    #     print('n=', n)
    #     print('inflow=', inflow)
    #     print('outflow=', outflow)
    #     print('sum(inflow)=', sum(inflow))
    #     print('sum(outflow)=', sum(outflow))
    # print('sum(recombined)=', sum(recombined))
    # print('recomb_helper=', recomb_helper)
    # print('n=', n)
    # plt.plot(e, n)
    # plt.xlabel('$E/\Delta$')
    # plt.ylabel('QP Density')
    # plt.yscale('log')
    # plt.show()
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
    idx = (np.abs(e - e_inj)).argmin()

    Vol_Al = 2.0  # um^3
    t = tau0 * dt * 1e-9  # in units of seconds
    N_qp = 2 * t * Gamma_p  # time * parity rate = QPs generated, 2 means a broken pair generates 2 QPs
    # print('N_qp=', N_qp)
    n_cp = 4.0 * 1e6  # Cooper pair density um^-3
    n_qp = N_qp / (Vol_Al)
    n_inj = (Delta / 2.0) * (n_qp / n_cp)
    delta_inj[idx] = n_inj  # number of QPs, not the reduced x_QP, this only linear scales the result

    return delta_inj


Gamma_p = 2*5.6e3 # at 250 GHz, one photon generates two QPs at half energy
n_inj = inj(ne, e_inj, Gamma_p, tau0, dt)
# print('n_inj=', n_inj)
### n_inj = [0,0,0,...,8.8e-10,...,0]

### evolve in time
n_before = n
i_converge = 0
for i in range(nt):
    # plt.plot(e, n)
    # plt.xlabel('$E/\Delta$')
    # plt.ylabel('QP Density')
    # plt.yscale('log')
    # plt.show()
    # if i % 1000 == 0:
    #     print(i/1000)
    # if i % 1000 == 0:

    # no spatial consideration here
    inj = n_inj
    n = scatter(n, inj, dt, i)
    # ChiSquare = Chi_Squaire(n, n_before)
    myC = myConverge(n, n_before)
    if i % 10000 == 0:
        print(i / 10000)
        # print('ChiSquare=', ChiSquare)
    if myC:
        print('Converge once')
        print('i_converge=', i_converge)
        i_converge = i_converge + 1
        if i_converge >= 5: break  # convergence criterion satisfied 5 times in a row
    else:
        i_converge = 0
    n_before = n

# print('sum(n[0, :])=', sum(n[0, :]))
print('x_qp=', 2 * sum(n))
ni = n

QP_Data = []
for i in range(len(ni)):
    d = [e[i], ni[i]]
    QP_Data.append(d)
np.savetxt('QPEnergySpectrumWithTrapping.txt', QP_Data)

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
