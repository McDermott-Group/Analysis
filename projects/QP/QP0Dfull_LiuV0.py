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
from __future__ import division
# from pylab import plot, xlabel, ylabel, show
import matplotlib.pyplot as plt
import random
import time
import numpy as np
from scipy.constants import *
from scipy.linalg import expm
from scipy.sparse import diags
import cmath


def myConverge(n_curr, n_before):
    myC = 1
    for i in range(len(n_curr) - 4):
        res = np.abs(n_curr[i] - n_before[i]) / n_curr[i]
        if res > 1.0e-10:  # or 1e-9
            # print('res=', res)
            myC = 0
            break
    return myC

Vol_Al = 6.4  # um^3
n_cp = 4.0 * 1e6

Delta = 1  # energy gap
# L = 2  # length in microns, assume in 2 um, the n_QP is uniform distribute
tau0 = 438  # e-p coupling time, units of ns, use Kaplan's value
# Dn = 6e-2  # units of um^2/ns -- use 6 here for aluminum

nx = 1  # points in space grid, let's say it is 1 D now
ne_log = 50  # points in log energy grid,
# 1000 enough for my purpose, 100 good for debug purpose
# finer ne also means finer time step

time_tot = 10 * 1e-3  # milli sec. 1 msec should be enough, from 1 to 10 msec,
# the change is 0.001% level, 100 msec is more than enough without trapping
s = 1e-3  # sets time step; units of tau0 -- try 1e-4, 5e-2 good for E=2Delta injection
dt = s * tau0  # time step, units of sec
Gamma = 1e-3  # Dynes broadening

nt = int(time_tot *1e9/ dt)  # time steps, convert things to nanoseconds
print('nt=', nt)
# ne = 100
# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
emax = 10.0  # units of Delta
# de = (emax - 1.0) / ne  # energy increment in units of Delta = 1
# e = 1.0 + de * np.linspace(0, ne - 1, ne)
e_phonon_inj = 5

### log grid
# e_log = np.logspace(0, 4, num=ne_log, base=2)
# idx_e_upper = (np.abs(e_log - e_inj)).argmin()
# e = e_log[:idx_e_upper + 5]  # useful energy grid

### uniform grid
e_uni = np.linspace(1, 10, 400)
e_uni = np.linspace(1, 10, 400)
idx_e_upper = (np.abs(e_uni - e_phonon_inj)).argmin()
e = e_uni[:idx_e_upper + 5]  # useful energy grid


ne = len(e)  # points in energy grid

### broadening is not the issue
broadened = complex(0, -Gamma)
Dynes = e + broadened  # to avoid divergence at the gap
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real

n = np.zeros(ne)  # n[i] is the normalized qp density per energy bin with units of energy

### initial data
if 0:  # use or not use initial data
    dat = np.loadtxt('QPEnergySpectrum_250GHzPhoton.txt')
    e = dat[:, 0]
    n = dat[:, 1]  # use earlier data as initial

G = np.zeros((ne, ne))  # scattering matrix
Gr0 = np.zeros((ne, ne))  # recombination matrix
Gp = np.zeros((ne, ne))  # phonon matrix

inflow = np.zeros(ne)
outflow = np.zeros(ne)
phonons_s = np.zeros(ne)

# generate the scattering and recombination matrices

### kTc = \Delta/1.76
for i in range(ne):
    for j in range(ne):
        # this is scattering from j to i
        phonon_factor = np.heaviside(e[j] - e[i], 0.5)
        ### There is a difference between the 0 Temp approximation and the temperature dependent one
        # if i == j:
        #     phonon_factor = 0   # to avoid divergence. there is ei-ej in other parts anyway
        # else:
        #     phonon_factor = (1.0/np.abs(np.exp(-(e[i] - e[j])/0.01)-1))

        # if e[j] < e[i]:
        #     phonon_factor = 0
        # else:
        #     phonon_factor = 1

        # print ('phonon_factor=', phonon_factor)

        # print('PhononTemp=', np.heaviside(e[j] - e[i], 0.5))
        # each time step = dt/tau0 = s
        G[i, j] = s * (1.76 ** 3.0) * ((e[i] - e[j]) ** 2.0) * (1.0 - Delta ** 2.0 / (e[i] * e[j])) \
                  * phonon_factor * rho[j]
        # if e[i] - e[j] > 0:
        #     G[i, j] = s * (1.76 ** 3.0) * ((e[i] - e[j]) ** 2.0) * (1.0 - Delta ** 2.0 / (e[i] * e[j])) \
        #               * rho[j]

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
column = np.ones(ne)
column.shape = (ne, 1)
outflow_helper = np.matmul(GT, column)


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

    # n = n + inflow - outflow - recombined-trapped
    n = n + inflow - outflow - recombined
    # n = n + inflow - outflow
    n = n + inj

    return n


# injection
def inj_delta(ne, e_inj, Gamma_qp, dt):
    """
    delta injection. at energy e_inj, with rate Gamma_qp and time step dt
    :param ne:  energy bins
    :param e_inj: delta injection energy, e.g. injection at 3*Delta
    :param Gamma_p: the measured parity rate in units of seconds, how many QPs are injected per second
    :return: inj, the injection array in units of Lenander paper, the same unit of ni
    """
    delta_inj = np.zeros(ne)
    delta_inj.shape = (ne, 1)
    idx = (np.abs(e - e_inj)).argmin()  # find the index of the delta injection energy

    N_qp = dt * Gamma_qp/1e9  # number of QP injected converted to each time step, dt in units of ns
    n_qp = N_qp/Vol_Al
    n_qp = n_qp * Delta / (2 * n_cp)
    # print('t=', t)
    # print('Gamma=', Gamma_qp/(Vol_Al*n_cp*Delta))
    print('n_qp=', n_qp)
    delta_inj[idx] = n_qp  # number of QPs, not the reduced x_QP, this only linear scales the result

    return delta_inj

def inj_dis(ne, e_p, Gamma_p, dt):
    """

    :param ne: points in energy grid
    :param e_p: inject phonon energy
    :param Gamma_p: injection rate units of Hz
    :param dt: time step in units of ns
    :return: distribution injection, same units as n
    """
    Np2q = np.zeros(ne)  # phonon to quasiparticle distribution, phonon number

    for p_i in range(ne):
        e_q = e[p_i]
        if 1 < e[p_i] < e_p -1:
            Np2q[p_i] = (e_q/np.sqrt(e_q**2-1)) * ((e_p-e_q)/np.sqrt((e_p-e_q)**2-1)) \
                        * (1 + 1/(e_q*(e_p-e_q)))

    p2qEnergyTot=sum(np.multiply(e, Np2q))
    normarlizeFactor = e_p/p2qEnergyTot # e_p is the total power injected at each step

    # print('p2qEnergyTot=', p2qEnergyTot)
    # print('normarlizeFactor=', normarlizeFactor)

    Np2q = Np2q * normarlizeFactor    # sum of p2q = e_p
    time_step = dt*Gamma_p/1e9
    Np2q = Np2q * time_step   # number

    np2q = Np2q / Vol_Al    # density per um^3
    np2q = np2q * Delta / (2 * n_cp)

    np2q.shape = (ne, 1)
    ### normalize energy

    # plt.plot(e, p2q)
    # plt.show()

    return np2q


# Gamma_qp = 1e9 #
Gamma_phonon = 1.0e9  #
Gamma_qp = Gamma_phonon  #
n_inj_delta = inj_delta(ne, e_phonon_inj, Gamma_qp, dt)
n_inj_dis = inj_dis(ne, e_phonon_inj, Gamma_qp, dt)
# plt.plot(e, n_inj_delta, label='Delta')
# plt.plot(e, n_inj_dis, label='Distributed')
# plt.xlabel('Energy/Delta')
# plt.ylabel('ni')
# plt.yscale('log')
# plt.legend()
# plt.show()


### evolve in time
n_before = n
i_converge = 0
convert_time_step = 0
n_inj = n_inj_dis
# n_inj = n_inj_delta
for i in range(nt):
    # for i in range(0):
    # no spatial consideration here
    if i % 10000 == 0:
        # print(i / 10000)
        print(int((i / nt)*100), '%')
    if i == 0:  # inject at the beginning
        inj = n_inj
    else:
        inj = n_inj * 1
    n = scatter(n, inj, dt, i)
    # myC = myConverge(n, n_before)
    # if myC:
    #     print('Converge once')
    #     print('i_converge=', i_converge)
    #     i_converge = i_converge + 1
    #     if i_converge >= 5:
    #         print('convert')
    #         convert_time_step = i
    #         break  # convergence criterion satisfied 5 times in a row
    # else:
    #     i_converge = 0
    n_before = n

print('x_qp=', 2*sum(n) / Delta)

QP_Data = []
energy_tot = 0.0
for i in range(len(n)):
    d = [e[i], n[i]]
    QP_Data.append(d)
    energy_tot = energy_tot + e[i] * n[i]

e_res = 0
e_inj = 0
for i in range(len(e)):
    e_res = e_res + (2*n_cp*Vol_Al/Delta)*e[i]*n[i]
    e_inj = e_inj + (2*n_cp*Vol_Al/Delta)*e[i]*n_inj[i]

e_inj = e_inj * nt
e_inj_direct_cal = e_phonon_inj * nt* dt * Gamma_phonon / 1e9

# e_res = np.multiply(e, n)
# e_inj = np.multiply(e, n_inj)

print('energy_res=', e_res)
print('injection energy from array=', e_inj)
print('injection energy direct calculation=', e_inj_direct_cal)
print('energy_res/energy_inj=', e_res / e_inj)
print('2Delta/e_inj=', 2/e[idx_e_upper])
print('convert_time_step=', convert_time_step)
np.savetxt('QP0D.txt', QP_Data)

plt.plot(e, n)
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
