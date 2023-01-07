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


Delta = 1 #energy gap
L = 100 #length in microns
tau0 = 400 #e-p coupling time, units of ns
Dn = 6e-2#units of um^2/ns -- use 6 here for aluminum


nx = 100 #points in space grid
ne = 300 #points in energy grid
nt = 100 #time steps

s = 1e-3 #sets time step; units of tau0 -- try 1e-4
dt = s*tau0
Gamma = 1e-3 # Dynes broadening

# spatial stuff
x = np.linspace(0,L,nx)
dx = x[1]-x[0]

# energy stuff
## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
emax = 2.0 #units of Delta
de = (emax-1.0)/ne #energy increment in units of Delta = 1
e = 1.0 + de*np.linspace(0,ne-1,ne)
broadened = complex(0,-Gamma)
Dynes = e + broadened
# rho = e/np.sqrt(e**2-Delta**2)
rho = Dynes/np.sqrt(Dynes**2-Delta**2)
rho = rho.real


n = np.zeros((nx,ne))

# print(e)

G = np.zeros((ne,ne))
Gr0 = np.zeros((ne,ne))
Gp = np.zeros((ne,ne))

inflow = np.zeros(ne)
outflow = np.zeros(ne)
phonons_s = np.zeros(ne)



#generate the scattering and recombination matrices
for i in range(ne):
    for j in range(ne):
        G[i,j]=s*(1.76**3)*((e[i]-e[j])**2)*np.heaviside(e[j]-e[i],0.5)*rho[i]*(1-Delta**2/e[i]/e[j])
        #this is scattering from j to i
        Gr0[i,j]=s*(1.76**3)*((e[i]+e[j])**2)*(1+Delta**2/e[i]/e[j])
        
#generate the scattering matrix for phonons
for i in range(ne):
    for j in range(ne):
        if i<=j:
            Gp[i,j] = G[j-i,j]


logG=np.log(G)
GT = G.transpose()
logGT=np.log(GT)
#print(GT)
column=np.ones(ne)
column.shape=(ne,1)
outflow_helper = np.matmul(GT,column)
# print(outflow_helper)



#generate the diffusion matrix
Dmat = diags([1, -2, 1], [-1, 0, 1], shape=(nx,nx)).toarray()
Dmat[0,0] = -1
Dmat[nx-1,nx-1] = -1 #hard boundary conditions

def diffuse(n,D,dx,dt):
    dn = D*dt/(dx**2)*np.matmul(Dmat,n) 
    n = n + dn
    return n

def scatter(n,inj,dt): 
    n.shape = (len(n),1)
    inflow = np.matmul(G,n) # this is scattering into the state
    #print(inflow)
    outflow = outflow_helper*n # elementwise multiplication 
    recomb_helper = np.kron(n,np.ones(ne))
    Gr = 2*Gr0*recomb_helper
    recombined = np.matmul(Gr,n)
    n = n + inflow - outflow - recombined
    n = n + inj
    n.shape = (1,len(n))
    return n
    
    
# injection
delta_inj = np.zeros(ne)
delta_inj.shape = (ne,1)
delta_inj[ne-50] = 1e-5

# used in case there is no injection
noinj = np.zeros(ne)
noinj.shape = (ne,1)


#set up the figure
# fig, ax = plt.subplots()
# distribution = ax.imshow(n, cmap='bwr')
# ax.set_title("QP Distribution")
# mngr = plt.get_current_fig_manager()
# mngr.window.setGeometry(750,100,640, 545)
# plt.ion()

# label_list = ["1.2","1.4","1.6","1.8"]
# ax.set_xticks([60,120,180,240])
# ax.set_xticklabels(label_list)
# plt.xlabel(r"Energy ($\Delta$)")
# plt.ylabel(r"Position ($\mu$m)")

#evolve in time
for i in range(nt):
    
    for j in range(nx):
        if j==nx//2: #inject QPs into one energy bin in center of the wire
            inj = delta_inj
        else:
            inj = noinj
        n[j,:] = scatter(n[j,:],inj,dt)
        
    for k in range(ne):
        D = Dn*np.sqrt(1-(Delta/e[k])**2)
        n[:,k] = diffuse(n[:,k],D,dx,dt)
    
    # if i%100 == 0:
    #     ax.imshow(np.log(n), cmap='bwr') #np.log(n)
    #     plt.pause(0.1)   

    

e_phon = e-1.0
e_recomb = 2.0 + np.linspace(0,2*ne-2,2*ne-1)

# look at spectrum of scattered phonons
phonons_s = np.matmul(Gp,n[nx//2,:])
phonons_r = np.zeros(2*ne-1)

#look at spectrum of recombination phonons
for k in range(2*ne-1):
    # print(k)
    # for l in range(k//2+1): # stop is k
    #     if k==300:
    #         print(l) 
    #     phonons_r[k] = phonons_r[k] + Gr0[l,k-l]*n[nx//2,l]*n[nx//2,k-l]
    l = 0
    while (l<=k//2) and (k//2+l<ne-1): #need to fix this -- treat even and odd cases separately
    # think this is not quite right -- indexing probably screwed up slightly at top of energy range
        if k%2 == 0: #k even, sum from main diagonal
            phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l]*n[nx//2,k//2-l]*n[nx//2,k//2+l]
        else:
            phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l+1]*n[nx//2,k//2-l]*n[nx//2,k//2+l+1]
        l+=1
    # else:
    #     for l in range(min(ne-1,k)-k//2+1):
    #         phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l+1]*n[nx//2,k//2-l]*n[nx//2,k//2+l+1]

e_phon = e-1


e_recomb = 2 + 2*(emax-1)*np.linspace(0,1,2*ne-1)

plot(e_phon,np.log(phonons_s),'b',e_recomb,np.log(phonons_r),'r',e,np.log(n[nx//2,:]),'g')
# plot(e_recomb,phonons_r,'r')
plot()
show()

####
#### here down for debugging
####


# ne=10
# Delta = 1

# # energy stuff
# ## IF NEED HIGHER RESOLUTION AT GAP, LOG SPACING FROM DELTA
# emax = 2 #units of Delta
# de = (emax-1)/ne #energy increment in units of Delta = 1
# e = 1 + de*np.linspace(0,ne-1,ne)

# Gamma = 1e-3 # Dynes broadening
# broadened = complex(0,-Gamma)
# Dynes = e + broadened
# # rho = e/np.sqrt(e**2-Delta**2)
# rho = Dynes/np.sqrt(Dynes**2-Delta**2)
# rho = rho.real

# s = 1e-4

# n = np.zeros(ne)
# n[1] = 0.01
# # n[1] = 0.03
# # n[2] = 0.01
# n[8] = 0.2

# n = n + 1e-6

# e_recomb = 2 + np.linspace(0,2*ne-2,2*ne-1)

# G = np.zeros((ne,ne))
# Gr0 = np.zeros((ne,ne))
# Gp = np.zeros((ne,ne))

# #generate the scattering and recombination matrices
# for i in range(ne):
#     for j in range(ne):
#         G[i,j]=s*(1.76**3)*((e[i]-e[j])**2)*np.heaviside(e[j]-e[i],0.5)*rho[i]*(1-Delta**2/e[i]/e[j])
#         #this is scattering from j to i
#         Gr0[i,j]=s*(1.76**3)*((e[i]+e[j])**2)*(1+Delta**2/e[i]/e[j])
        
# logGr0=np.log(Gr0)

# # #set up the figure
# # fig, ax = plt.subplots()
# # QPs = ax.imshow(logGr0, cmap='bwr')
# # ax.set_title("QP Recombination Matrix")
# # mngr = plt.get_current_fig_manager()
# # mngr.window.setGeometry(750,100,640, 545)

#     #     plt.pause(0.1)   
        
# #generate the scattering matrix for phonons
# for i in range(ne):
#     for j in range(ne):
#         if i<=j:
#             Gp[i,j] = G[j-i,j]

# #look at spectrum of scattered phonons
# phonons_s = np.matmul(Gp,n)
# phonons_r = np.zeros(2*ne-1)

# #look at spectrum of recombination phonons
# # phonons_r[0] = Gr0[0,0]*n[0]**2
# for k in range(2*ne-1):
#     # print(k)
#     # for l in range(k//2+1): # stop is k
#     #     if k==300:
#     #         print(l) 
#     l = 0
#     # misses the 00 element since l=k here
#     while (l<=k//2) and (k//2+l<ne-1): #need to fix this -- treat even and odd cases separately
#         if k%2 == 0: #k even, sum from main diagonal
#             phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l]*n[k//2-l]*n[k//2+l]
#         else:
#             phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l+1]*n[k//2-l]*n[k//2+l+1]
#         l+=1
#     print(phonons_r[k])
#     print(l)
#     # else:
#     #     for l in range(min(ne-1,k)-k//2+1):
#     #         phonons_r[k] = phonons_r[k] + Gr0[k//2-l,k//2+l+1]*n[nx//2,k//2-l]*n[nx//2,k//2+l+1]

# e_phon = e-1


# e_recomb = 2 + 2*(emax-1)*np.linspace(0,1,2*ne-1)

# # plot(e_phon,np.log(phonons_s),'b',e_recomb,np.log(phonons_r),'r',e,np.log(n),'g')
# # plot(e,n,'r')
# plot(np.log(phonons_r),'r')
# plot()
# show()




        

   
    
