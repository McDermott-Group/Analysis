import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, UpAndParity, getPhotonRate
import numpy as np
import matplotlib.font_manager as font_manager
from scipy import interpolate
from scipy.integrate import quad

if 1:   # for BB induced uprate
    """
    Import CST files and Junction parameters
    """
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff_circ = 75 * 1e-21  # Commonly used (50-100)

    JQ4 = [15 * 1e3, 19.9 * 1e-9, 0, 184.4 * 122.5 * 2, "Receiver"]  #

    ### up to 1000 GHz

    fileQ4 = "Q4_full-chip.txt"
    Q4_Antenna = AntennaCoupling()
    Q4_Antenna.import_data(fileQ4, JQ4, C_eff=C_eff_circ)
    f_Q4 = Q4_Antenna.Antenna["f"]
    eQ4 = Q4_Antenna.Antenna["e_c"]
    ecQ4 = Q4_Antenna.Antenna["e_c_dB"]
    AreaQ4 = Q4_Antenna.Receiver["Area"]

if 1:
    """
    Import measurement data starts
    """
    if 1:
        """AugCircRadiator Data"""
        f_DAC_Aug = 4.604
        Q4_PSD_file_Aug = "2021AugCircRadiator_Q4_PSD_Data.txt"
        Q4_P1_file_Aug = "2021AugCircRadiator_Q4_P1_Data.txt"
        Q4_Up_file_Aug = "2021AugCircRadiator_Q4_Up_Data.txt"
        Q4_Up_1D_283 = "2021Sep01CircRadiator_Q4_Up_1D_Data_283Bias.txt"
        Q4_Up_1D_322 = "2021Sep01CircRadiator_Q4_Up_1D_Data_322Bias.txt"
        # Up rate for 283 GammaUp=543.6357980508438
        # Up rate for 322 GammaUp=624.5424275627762
        UpRate322 = 624.5424275627762

        Q4_Up_1D = np.loadtxt(Q4_Up_1D_322, skiprows=20)
        Up_1D_Time = Q4_Up_1D[:, 0]
        Up_1D_P1 = Q4_Up_1D[:, 1]
        Up_1D_Std = Q4_Up_1D[:, 2]
        Up_1D_P1Fit = Q4_Up_1D[:, 3]

        S_param_file = "S_minusplus.txt"
        #S_param_file = "S_minusplus_1THz.txt"

        """Working on the theory"""
        UpParity = UpAndParity()
        UpParity.import_data(Q4_PSD_file_Aug, Q4_Up_file_Aug, S_param_file)
        # Up_Base = UpParity.getBaseUpYale(antenna=Q4_Antenna)
        # print('Up_Base=', Up_Base)
        BB_UpRate = 83.7  # Hz, calculated from effective BB temp and Yale's theory

        parity_freq_data = np.array(UpParity.parity_freq_data)
        parity_data = np.array(UpParity.parity_data)
        up_freq_data = np.array(UpParity.up_freq_data)
        up_data = np.array(UpParity.up_data)
        UpRateYale = np.array(UpParity.UpRateYale)

        f_interest = np.array(UpParity.f_interest)
        parity = np.array(UpParity.parity_interest)
        up = np.array(UpParity.up_interest)
        parity_PAT = np.array(UpParity.parity_PAT)
        parity_QPD = np.array(UpParity.parity_QPD)
        up_PAT = np.array(UpParity.up_PAT)
        up_QPD = np.array(UpParity.up_QPD)

        ### plot start
        label_font = 22
        tick_font = 20
        legend_font = 14

        fig, axs = plt.subplots(2, figsize=(10, 8),
                                gridspec_kw={'height_ratios': [3, 4],
                                             'hspace': 0.65}
                                )

        axs[0].errorbar(Up_1D_Time, Up_1D_P1*100, yerr=100*Up_1D_Std/np.sqrt(20), fmt='o', c='k')
        axs[0].plot(Up_1D_Time, Up_1D_P1Fit*100, c='k')
        axs[0].set_xlim([3.8, 15.2])
        axs[0].set_ylim([1.15, 2.01])
        axs[0].set_yticks([1.2, 1.4, 1.6, 1.8, 2.0])

        axs[0].set_ylabel('P$_{1}$ (%)', fontsize=label_font, family='Arial')
        axs[0].tick_params(labelsize=tick_font)
        axs[0].tick_params(axis="x", direction="in", which='both')
        axs[0].tick_params(axis="y", direction="in", which='both')
        axs[0].tick_params(axis="x", width=1, length=4, which='both')
        axs[0].tick_params(axis="y", width=1, length=3, which='minor')
        axs[0].tick_params(axis="y", width=1, length=6, which='major')

        f_l = 0
        f_r = 620
        loc = 2

        start = 26

        # axs[1].plot(parity_freq_data, (parity_data-600)*0.22, 'r', linewidth=3, label='         ')
        # axs[1].plot(up_freq_data, up_data-300, 'b', linewidth=3, label='         ')
        # axs[1].plot(parity_freq_data[start:], UpRateYale[start:], 'b--', linewidth=3, label='         ')

        axs[1].plot(parity_freq_data, parity_data, 'r', linewidth=3, label='         ')
        axs[1].plot(up_freq_data, up_data, 'b', linewidth=3, label='         ')
        axs[1].plot(parity_freq_data[start:], UpRateYale[start:], 'b--', linewidth=3, label='         ')

        # axs[1].plot(f_interest, parity, 'r', linewidth=3, label='             ')
        # axs[1].plot(f_interest, up, 'b', linewidth=3, label='             ')
        # axs[1].plot(f_interest, up_PAT, 'b--', linewidth=3, label='             ')
        axs[1].set_ylim([6.5e1, 1.2e4])
        axs[1].set_yscale('log')
        axs[1].set_ylabel('Rate ($\mathrm{s}^{-1}$)', fontsize=label_font, family='Arial')
        axs[1].tick_params(labelsize=tick_font)
        # axs[1].legend(loc=loc, prop={'size': 13}, frameon=False)
        axs[1].set_xlim([f_l, f_r])
        # axs[1].set_xlabel('Frequency (GHz)', fontsize=label_font, family='Arial')

        axs[1].tick_params(axis="x", direction="in", which='both')
        axs[1].tick_params(axis="y", direction="in", which='both')
        axs[1].tick_params(axis="x", width=1, length=4, which='both')
        axs[1].tick_params(axis="y", width=1, length=3, which='minor')
        axs[1].tick_params(axis="y", width=1, length=6, which='major')

        axs[1].axvspan(368, 620, facecolor='0.2', alpha=0.1)


        def fj_to_delta(f):
            return f / 92


        def delta_to_fj(delta):
            return 92 * delta


        secax = axs[1].secondary_xaxis('top', functions=(fj_to_delta, delta_to_fj))
        secax.tick_params(labelsize=tick_font)
        secax.set_xticks(np.arange(0, 10, 1))
        secax.set_xlabel("Transmitter Voltage Bias ($\Delta/e$)", fontsize=label_font, labelpad=10)
        secax.tick_params(axis="x", direction="in", which='both')
        secax.tick_params(axis="y", direction="in", which='both')

        secax.tick_params(axis="x", width=1, length=6, which='both')
        secax.tick_params(axis="y", width=1, length=3, which='minor')
        secax.tick_params(axis="y", width=1, length=6, which='major')

        fig.align_ylabels(axs)
        plt.tight_layout()
        # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        # plt.savefig(path + '\ParityUpCorrelation.pdf', format='pdf', bbox_inches='tight', dpi=1200)


        plt.savefig('ParityUpCorrelation.pdf', format='pdf', bbox_inches='tight', dpi=1200)

        plt.show()

