from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os


### parameters to be tuned
e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
C_eff = 75 * 1e-21  # Commonly used (50-100) 75, 88
JQ2 = [16.6 * 1e3 * 100, None, 0, 360 * 150, "Receiver"]  #

fileQ2 = "xmon_full-chip_Q2.txt"

Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2, C_eff=C_eff)
f_Q2 = Q2.Antenna["f"]
eQ2 = Q2.Antenna["e_c"]
Z_Im = Q2.Antenna["Z_Im"]
Z_Re = Q2.Antenna["Z_Re"]
Z_j_star = np.conj(Q2.Junction["Z_j"])
Area = Q2.Receiver["Area"]

f_Q2 = f_Q2/1e9

label_font = 20
tick_font = 20
legend_font = 20
ld2 = 5

Z_j_star_real = Z_j_star.real
Z_j_star_imag = Z_j_star.imag

l_i = 20
r_i = 1000
plt.rcParams.update({'font.family': 'Arial'})

phfont = {'fontname': 'Times New Roman'}
# plt.rcParams['font.size'] = 12
# plt.rcParams.update({
#   "text.usetex": True,
#   "font.family": "Times New Roman"
# })
if 1:
    fig, axs = plt.subplots(3, sharex='col', figsize=(7.5, 8),
                            gridspec_kw={'hspace': 0.15})

    legend_font = font_manager.FontProperties(family='Times New Roman',
                                       weight='bold',
                                       style='normal', size=legend_font)

    axs[0].plot(f_Q2, Z_Re, linewidth=ld2, color="blue", label='$Z_{\\rm rad}$')
    axs[0].plot(f_Q2, Z_j_star_real, linewidth=ld2, color="red", label='$Z_{j}^*$')
    axs[0].set_xlim([l_i, r_i])
    axs[0].set_ylabel('Re[Z] ($\Omega$)', fontsize=label_font, **phfont)
    axs[0].set_ylim([0, 350])
    axs[0].set_yticks([0, 100, 200, 300])
    axs[0].tick_params(labelsize=tick_font)
    # axs[0].legend(loc=2, prop=legend_font, frameon=False)

    axs[0].tick_params(axis="x", direction="in", which='both')
    axs[0].tick_params(axis="y", direction="in", which='both')
    axs[0].tick_params(axis="x", width=1, length=6, which='both')
    axs[0].tick_params(axis="y", width=1, length=6, which='both')
    # plt.rcParams.update({'font.family': 'Arial'})

    axs[1].plot(f_Q2, Z_Im, linewidth=ld2, color="blue", label='$Z_{\\rm rad}$')
    axs[1].plot(f_Q2, Z_j_star_imag, linewidth=ld2, color="red", label='$Z_{j}^*$')
    axs[1].set_xlim([l_i, r_i])
    axs[1].set_ylim([-150, 200])
    axs[1].set_yticks([-100, 0, 100, 200])
    axs[1].tick_params(labelsize=tick_font)
    axs[1].set_ylabel('Im[Z] ($\Omega$)', fontsize=label_font, **phfont)
    axs[1].tick_params(axis="x", direction="in", which='both')
    axs[1].tick_params(axis="y", direction="in", which='both')
    axs[1].tick_params(axis="x", width=1, length=6, which='both')
    axs[1].tick_params(axis="y", width=1, length=6, which='both')
    # axs[1].legend(loc=2, prop={'size': legend_font}, frameon=False)

    axs[2].plot(f_Q2, eQ2, linewidth=ld2, color="black")
    axs[2].set_xlim([l_i, r_i])
    axs[2].set_ylim([2e-4, 1e-1])
    axs[2].tick_params(labelsize=tick_font)
    axs[2].set_yscale('log')
    axs[2].set_ylabel('$e_{\mathrm{c}}$', fontsize=label_font, **phfont)
    axs[2].set_xlabel('Frequency (GHz)', fontsize=label_font)
    axs[2].set_xticks(np.arange(100, 1000, 200))
    axs[2].tick_params(axis="x", direction="in", which='both')
    axs[2].tick_params(axis="y", direction="in", which='both')
    axs[2].tick_params(axis="x", width=1, length=6, which='both')
    axs[2].tick_params(axis="y", width=1, length=6, which='major')
    axs[2].tick_params(axis="y", width=1, length=3, which='minor')

    fig.align_ylabels(axs[:])
    plt.tight_layout()

    # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    # plt.savefig(path+'\AntennaImpedance.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

if 0:   # zoom in
    fig = plt.figure(figsize=(5, 3))
    plt.plot(f_Q2, Z_Re, linewidth=ld2, color="blue", label='$Z_{\\rm rad}$')
    plt.plot(f_Q2, Z_j_star_real, linewidth=ld2+2, color="red", label='$Z_{j}^*$')
    plt.xlim([l_i, r_i])
    plt.xlim([150, 250])
    plt.ylim([0, 7.1])

    plt.tick_params(axis="x", direction="in", which='both')
    plt.tick_params(axis="y", direction="in", which='both')
    plt.tick_params(axis="x", width=2, length=10, which='both')
    plt.tick_params(axis="y", width=2, length=10, which='both')

    plt.xticks([160, 200, 240])
    plt.yticks([1, 3, 5, 7])
    # plt.xticklabels(f_range)

    plt.tick_params(labelsize=36)


    path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    plt.savefig(path+'\AntennaImpedance_Zoomin.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()
