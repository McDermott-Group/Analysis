#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:00:25 2021

@author: robertmcdermott

Modified:
2022Jan
Chuanhong Liu

For recombination photon emission

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
ne = 300  # points in energy grid
nt = 10  # time steps

s = 1e-3  # sets time step; units of tau0 -- try 1e-4
dt = s * tau0
Gamma = 1e-3  # Dynes broadening


# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
emax = 2.0  # units of Delta, this needs to be modified for different bias voltage
de = (emax - 1.0) / ne  # energy increment in units of Delta = 1
e = 1.0 + de * np.linspace(0, ne - 1, ne)
broadened = complex(0, -Gamma)
Dynes = e + broadened
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real

n = np.zeros((nx, ne))

G = np.zeros((ne, ne))
Gr0 = np.zeros((ne, ne))
Gp = np.zeros((ne, ne))

inflow = np.zeros(ne)
outflow = np.zeros(ne)
phonons_s = np.zeros(ne)

# generate the scattering and recombination matrices
for i in range(ne):
    for j in range(ne):
        G[i, j] = s * (1.76 ** 3) * ((e[i] - e[j]) ** 2) * np.heaviside(e[j] - e[i], 0.5) * rho[i] * (
                    1 - Delta ** 2 / e[i] / e[j])
        # this is scattering from j to i
        Gr0[i, j] = s * (1.76 ** 3) * ((e[i] + e[j]) ** 2) * (1 + Delta ** 2 / e[i] / e[j])

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
    n.shape = (len(n), 1)
    inflow = np.matmul(G, n)  # this is scattering into the state
    # print(inflow)
    outflow = outflow_helper * n  # elementwise multiplication
    recomb_helper = np.kron(n, np.ones(ne))
    Gr = 2 * Gr0 * recomb_helper
    recombined = np.matmul(Gr, n)
    n = n + inflow - outflow - recombined
    n = n + inj
    n.shape = (1, len(n))
    return n

# injection
delta_inj = np.zeros(ne)
delta_inj.shape = (ne, 1)
delta_inj[ne - 50] = 1e-5

# used in case there is no injection
noinj = np.zeros(ne)
noinj.shape = (ne, 1)

# evolve in time
for i in range(nt):
    # no spatial consideration here
    inj = delta_inj
    n[0, :] = scatter(n[0, :], inj, dt)

e_phon = e - 1.0
e_recomb = 2.0 + np.linspace(0, 2 * ne - 2, 2 * ne - 1)

# look at spectrum of scattered phonons
phonons_s = np.matmul(Gp, n[nx // 2, :])
phonons_r = np.zeros(2 * ne - 1)

# look at spectrum of recombination phonons
for k in range(2 * ne - 1):
    l = 0
    while (l <= k // 2) and (k // 2 + l < ne - 1):  # need to fix this -- treat even and odd cases separately
        # think this is not quite right -- indexing probably screwed up slightly at top of energy range
        if k % 2 == 0:  # k even, sum from main diagonal
            phonons_r[k] = phonons_r[k] + Gr0[k // 2 - l, k // 2 + l] * n[nx // 2, k // 2 - l] * n[nx // 2, k // 2 + l]
        else:
            phonons_r[k] = phonons_r[k] + Gr0[k // 2 - l, k // 2 + l + 1] * n[nx // 2, k // 2 - l] * n[
                nx // 2, k // 2 + l + 1]
        l += 1

e_phon = e - 1
e_recomb = 2 + 2 * (emax - 1) * np.linspace(0, 1, 2 * ne - 1)

# plot(e_phon, np.log(phonons_s), 'b', e_recomb, np.log(phonons_r), 'r')
plt.plot(e_phon, phonons_s, 'b', label='scatter')
plt.plot(e_recomb, phonons_r, 'r', label='recomb')
plt.yscale('log')
plt.legend()
plt.show()
