#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:00:25 2021

@author: robertmcdermott
"""

from pylab import plot, xlabel, ylabel, show
import matplotlib.pyplot as plt
import random
import time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
import cmath

Delta = 1  # energy gap
L = 100  # length in microns
tau0 = 400  # e-p coupling time, units of ns
Dn = 6e-2  # units of um^2/ns -- want to use 6 here for aluminum

nx = 100  # points in space grid
ne = 300  # points in energy grid
# nt = 400  # time steps
nt = 100  # time steps

s = 1e-3  # sets time step; units of tau0 -- try 1e-4
dt = s * tau0
Gamma = 1e-3  # Dynes broadening

# spatial stuff
x = np.linspace(0, L, nx)
dx = x[1] - x[0]

# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
de = 1 / ne  # energy increment in units of Delta = 1
e = 1 + de * np.linspace(0, ne - 1, ne)
broadened = complex(0, -Gamma)
Dynes = e + broadened
# rho = e/np.sqrt(e**2-Delta**2)
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real

n = np.zeros((nx, ne))

# print(e)

G = np.zeros((ne, ne))
Gr0 = np.zeros((ne, ne))

inflow = np.zeros(ne)
outflow = np.zeros(ne)

# generate the scattering and recombination matrices
for i in range(ne):
    for j in range(ne):
        G[i, j] = s * (1.76 ** 3) * ((e[i] - e[j]) ** 2) * np.heaviside(e[j] - e[i], 0.5) * rho[i] * (
                    1 - Delta ** 2 / e[i] / e[j])
        # this is scattering from j to i
        Gr0[i, j] = s * (1.76 ** 3) * ((e[i] + e[j]) ** 2) * (1 + Delta ** 2 / e[i] / e[j])

logG = np.log(G)
GT = G.transpose()
logGT = np.log(GT)
# print(GT)
column = np.ones(ne)
column.shape = (ne, 1)
outflow_helper = np.matmul(GT, column)
# print(outflow_helper)


# generate the diffusion matrix
Dmat = diags([1, -2, 1], [-1, 0, 1], shape=(nx, nx)).toarray()
Dmat[0, 0] = -1
Dmat[nx - 1, nx - 1] = -1


def diffuse(n, D, dx, dt):
    dn = D * dt / (dx ** 2) * np.matmul(Dmat, n)
    n = n + dn
    return n


def scatter(n, inj, dt):
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
delta_inj[ne - 10] = 1e-5

# used in case there is no injection
noinj = np.zeros(ne)
noinj.shape = (ne, 1)

# set up the figure
fig, ax = plt.subplots()
distribution = ax.imshow(n, cmap='bwr')
ax.set_title("QP Distribution")
# mngr = plt.get_current_fig_manager()
# mngr.window.setGeometry(750, 100, 640, 545)
plt.ion()

label_list = ["1.2", "1.4", "1.6", "1.8"]
ax.set_xticks([60, 120, 180, 240])
ax.set_xticklabels(label_list)
plt.xlabel(r"Energy ($\Delta$)")
plt.ylabel(r"Position ($\mu$m)")

# evolve in time
for i in range(nt):
    print('i=', i)

    for j in range(nx):
        if j == nx // 2:
            inj = delta_inj
        else:
            inj = noinj
        n[j, :] = scatter(n[j, :], inj, dt)

    for k in range(ne):
        D = Dn * np.sqrt(1 - (Delta / e[k]) ** 2)
        n[:, k] = diffuse(n[:, k], D, dx, dt)

    if i % 100 == 99:
        print('i=', i)
        print('plot movie')
        ax.imshow(np.log(n), cmap='bwr')  # np.log(n)
        plt.pause(0.1*100)

ax.imshow(np.log(n), cmap='bwr')  # np.log(n)
plt.show()









