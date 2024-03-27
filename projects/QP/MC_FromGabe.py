from __future__ import division

import numpy as np
import numpy.random as rand


import matplotlib.pyplot as plt
import scipy.integrate as scint

hbar = 6.5821e-16 #eV s
e = 1.602e-19 #Coulombs - charge of an electron
phi_0 = 2.07e-15 #Vs magnetic flux quantum
k = 8.617e-5 # eV/K Boltzmanns constant in eV/K
Delta_eV =190e-6 #eV superconducting gap energy of Al
Delta_K = Delta_eV/k #K superconducting gap energy of Al
# n_cp = 2.8e6 # 1/um^3 number density of cooper pairs
n_cp = 4e6 # 1/um^3 number density of cooper pairs


def rho(epsilon, Delta):
    """
    Parameters:
    epsilon - energy of state
    Delta - superconducting gap energy

    returns:
    rho from equation 2 of the paper
    """
    return epsilon / (np.sqrt(epsilon ** 2 - Delta ** 2))


def P_qp(E_p, e, Delta):
    """
    Parameters:
    E_p - initial phonon energy
    e - Energy of one of the quasiparticles
    Delta - superconducting gap energy

    returns:
    Distribution of the resulting QP pair
    """
    # dealing with divide by zero errors
    if (E_p - e) ** 2 - Delta ** 2 < 1e-9:
        return 0
    if (e) ** 2 - Delta ** 2 < 1e-9:
        return 0
    return rho(e, Delta) * rho((E_p - e), Delta) * (1 + (Delta ** 2 / (e * (E_p - e))))


def P_ph(E_q, E_new, Delta):
    """
    Parameters:
    E_q - initial quasiparticle energy
    E_new - new quasiparticle energy
    Delta - superconducting gap energy

    returns:
    Distribution of the resulting phonon from scattering
    """
    # dealing with divide by zero errors
    if (E_new) ** 2 - Delta ** 2 < 1e-9:
        return 0
    kTc = 1 / 1.76  # k T_c for an Al junction
    return ((E_q - E_new) ** 2) * (1 / (kTc ** 3)) * rho(E_new, Delta) * (1 - (Delta ** 2) / (E_q * E_new))


def MC_Energy_Branching(E_p_0, Delta, iters):
    """
    Parameters:
    E_p_0 - initial phonon energy
    Delta - superconducting gap energy
    iters - number of Monte Carlo Schemes to average over

    returns:
    average *number* of resulting quasiparticles from the energy cascade
    """
    i = 0
    nbins = 100  # Number of Riemann bins for integrating PDFs

    # list of phonon energies
    E_phonons = [E_p_0]
    # list of qp energies
    E_qp = []

    list_qp_averages = []
    while i < iters:
        # End condition: all phonons are no longer pair breaking, all quasiparticles cannot create pair-breaking phonons
        if all(val <= 2 * Delta for val in E_phonons) and all(val2 <= 3 * Delta for val2 in E_qp):
            # record the *number* of qps created this iteration
            list_qp_averages.append(len(E_qp))
            # Resetting the lists...
            E_phonons = [E_p_0]
            E_qp = []
            i += 1
            continue

        else:
            # Pair-breaking
            for i_ph in np.arange(len(E_phonons)):
                E_ph_0 = E_phonons[i_ph]
                # Check for adequate pair-breaking energy
                if E_ph_0 < 2 * Delta:
                    continue
                # The approach is to randomly sample a quasiparticle energy from pair-break PDF by integrating
                # into a numerical CDF curve, inverting that curve, and plug in a random number

                # create CDF lookup table
                ltable = np.zeros(nbins)
                # calculate bin size
                bsize = (E_ph_0 - 2 * Delta) / nbins
                # Create a rough CDF from the PDF using Riemann integration
                total = 0
                for j in range(nbins):
                    riemann_height = P_qp(E_ph_0, (Delta + (j + 1 / 2) / nbins * (E_ph_0 - 2 * Delta)), Delta)
                    # Sanity check
                    if riemann_height < 0:
                        print("Pair-breaking PDF is returning negative value!")
                    ltable[j] = total + bsize * riemann_height
                    total += bsize * riemann_height
                # normalize
                ltable = ltable / ltable[(nbins - 1)]
                # random number from 0 to 1
                u = rand.uniform()
                # Find what bin value this is closest to
                bestbin = (np.abs(ltable - u)).argmin()

                # map this random number to the energy of (one of) the new qp
                E_qp_new = Delta + (bestbin + 1 / 2) / nbins * (E_ph_0 - 2 * Delta)

                # We now have two quasiparticles
                E_qp.append(E_qp_new)
                E_qp.append(E_ph_0 - E_qp_new)

                # our phonon is also gone
                E_phonons[i_ph] = 0

            # Quasiparticle Scattering
            for i_qp in np.arange((len(E_qp))):
                E_qp_0 = E_qp[i_qp]
                # Check if QP has enough energy to create pairbreaking phonon
                if E_qp_0 < 3 * Delta:
                    continue

                # Same random sampling approach as above

                # create lookup table
                ltable = np.zeros(nbins)
                # calculate bin size
                bsize = (E_qp_0 - Delta) / nbins
                # rough cdf
                total = 0
                for j in range(nbins):
                    riemann_height = P_ph(E_qp_0, (Delta + (j + 1 / 2) / nbins * (E_qp_0 - Delta)), Delta)
                    if riemann_height < 0:
                        print("Scattering PDF is returning negative value!")

                    ltable[j] = total + bsize * riemann_height
                    total += bsize * riemann_height
                # normalize
                ltable = ltable / ltable[nbins - 1]

                # random number from 0 to 1
                u = rand.uniform()

                # Find what bin value this is closest to
                bestbin = (np.abs(ltable - u)).argmin()

                # map this random number to the energy of the new qp
                E_qp_new = Delta + (bestbin + 1 / 2) / nbins * (E_qp_0 - Delta)

                # We now have a new phonon
                E_phonons.append(E_qp_0 - E_qp_new)
                # Our quasiparticle has a lower energy now!
                E_qp[i_qp] = E_qp_new

    return np.average(list_qp_averages)

# print('Here')
# print(MC_Energy_Branching(10 * Delta_K, Delta_K, 500) / 10)

E_ph = np.linspace(0.01*Delta_K,14*Delta_K,10)
n_qp = np.zeros(len(E_ph))
for i in range(len(E_ph)):
    print(i)
    n_qp[i] = MC_Energy_Branching(E_ph[i],Delta_K,1000)/(E_ph[i]/Delta_K)
    print(n_qp)
fig = plt.figure(figsize=(12,8), dpi= 100, facecolor='w', edgecolor='k')
ax = plt.gca()
plt.scatter(E_ph/Delta_K,n_qp, c = 'k')
plt.axhline(y=0.57,color = 'r',label = r'57% Conversion')
plt.title(r' Phonon to Quasiparticle down conversion (Fig 4)')
plt.xlabel(r'Initial phonon Energy $(E_p/\Delta)$')
plt.ylabel(r' Quasiparticle fraction $n_{qp} / (E_p/\Delta)$')
plt.legend()
plt.show()