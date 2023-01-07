#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:00:25 2021

@author: robertmcdermott

Modified:
2022Jan
Chuanhong Liu

For recombination photon emission

Step:
1. Inject
2. Recombine
3. Diffuse

Ignore the scattering process, only energy at the gap,
"""

from pylab import plot, xlabel, ylabel, show
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags
import cmath

Delta = 1  # energy gap
L = 40  # length in microns
# L = 4  # length in microns
tau0 = 400  # e-p coupling time, units of ns
Dn = 6e-2  # units of um^2/ns -- want to use 6 here for aluminum

nx = 100  # points in space grid
ne = 1  # points in energy grid
nt = 1000# time steps

s = 1e-3  # sets time step; units of tau0 -- try 1e-4
dt = s * tau0
Gamma = 1e-3  # Dynes broadening

# spatial stuff
x = np.linspace(0, L, nx)
dx = x[1] - x[0]

# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
de = 1.0 / ne  # energy increment in units of Delta = 1
e = np.asarray([1.05])
broadened = complex(0, -Gamma)
Dynes = e + broadened
# rho = e/np.sqrt(e**2-Delta**2)
rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
rho = rho.real

# generate the diffusion matrix
Dmat = diags([1, -2, 1], [-1, 0, 1], shape=(nx, nx)).toarray()
if 1:   # QPs confinement, without this ,there will be leakage from the bonds to Nb
    Dmat[0, 0] = -1
    Dmat[nx - 1, nx - 1] = -1

def diffuse(n, D, dx, dt):
    dn = D * dt / (dx ** 2) * np.matmul(Dmat, n)
    n = n + dn
    return n

def recombine(n, inj, dt, dx, i, j):
    n.shape = (len(n), 1)
    Gamma_rec = 1.0/(tau0)   # recombination rate, use the same rate as e-phonon coupling rate
    # Gamma_rec = 1.0/(100)   # recombination rate
    # absolute number of QPs to density conversion
    Vol = dx * 1.6 * 0.1    # Volume per segment    um^3
    n_cp = 4.0e6    # cooper pair density um^(-3)
    N_cp = n_cp*Vol   # cooper pair number per unit length
    x_qp = n/(1.0*N_cp)   # normalized QP density per unit length
    recombined = Gamma_rec*dt*x_qp*x_qp   # elementwise multiply
    x_qp = x_qp - recombined
    n = x_qp * N_cp
    n_rec = recombined * N_cp
    n = n + inj
    n.shape = (1, len(n))
    return n, n_rec, x_qp

# injection
E = 1.602*10**(-19) # elementary charge
Rn = 33.0e3   # resistance Ohm
# Vb = 1.45e-3  # bias voltage Volt
Vb = 2.0e-3  # bias voltage Volt
Delta_0 = 190*E*1e-6    # energy gap
G_inj = 0.57*(Vb**2/Rn)*(1/Delta_0)
G_inj_dt = G_inj*dt*1e-9    # injection QP number per time step
# print('G_inj_dt=', G_inj_dt)

delta_inj = np.zeros(ne)
delta_inj.shape = (ne, 1)
delta_inj[0] = G_inj_dt

# used in case there is no injection
noinj = np.zeros(ne)
noinj.shape = (ne, 1)

# evolve in time
n = np.zeros((nx, ne))
x_qp = np.zeros((nx, ne))
n_rec = np.zeros((nx, ne))
for i in range(nt):
    for j in range(nx): # spatial step
        if j == nx // 2:    # only injection at the center
            inj = delta_inj
        else:
            inj = noinj
        n[j], n_rec[j], x_qp[j] \
            = recombine(n[j], inj, dt, dx, i, j)

    D = Dn * np.sqrt(1 - (Delta / e[0]) ** 2)
    n = diffuse(n, D, dx, dt)

    photon_emitted = np.sum(n_rec) / 2.0

    if i % 5000 == 0:
        print('i=', i)
        # print('n_rec=', n_rec)
        print('photon_emitted=', photon_emitted)
        # print('np.sum(n)-sum_temp=', np.sum(n)-sum_temp)
        # sum_temp = np.sum(n)

print('Vb=', Vb)
# plt.plot(x, n)
plt.plot(x, x_qp)
plt.xlabel('Al lead (um)')
plt.ylabel('x_qp')
plt.show()

# rec_emitted = np.array([[0.4, 8.81], [0.5, 16.1], [0.6, 25.9], [0.7, 38.0], [0.8, 52.6], [0.9, 69.5],
#                [1.0, 88.9], [1.1, 110], [1.2, 135], [1.3, 161], [1.4, 190], [1.5, 221], [1.6, 255],
#                [1.7, 290], [1.8, 329], [1.9, 369], [2.0, 412]])
#
# plt.plot(rec_emitted[:, 0], rec_emitted[:, 1], label='photon emitted')
# # plt.plot(rec_emitted[:, 0], 100*rec_emitted[:, 0]*rec_emitted[:, 0], label='square dependence')
# plt.plot(rec_emitted[:, 0], 100*rec_emitted[:, 0]**2.0, label='square dependence')
# plt.plot(rec_emitted[:, 0], 25*rec_emitted[:, 0]**4.0, label='quadratic dependence')
# plt.xlabel('Vb (mV)')
# plt.ylabel('Photon Emitted per unit time')
# plt.legend()
# plt.show()









