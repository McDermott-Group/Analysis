import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, UpAndParity, getPhotonRate
import numpy as np
import matplotlib.font_manager as font_manager
from scipy import interpolate
from scipy.integrate import quad

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

        """Working on the theory"""
        UpParity = UpAndParity()
        UpParity.import_data(Q4_PSD_file_Aug, Q4_Up_file_Aug, S_param_file)
        BB_UpRate = 52.7  # Hz, calculated from effective BB temp and Yale's theory

        parity_freq_data = np.array(UpParity.parity_freq_data)
        parity_data = np.array(UpParity.parity_data)
        up_freq_data = np.array(UpParity.up_freq_data)
        up_data = np.array(UpParity.up_data)

        f_interest = np.array(UpParity.f_interest)
        parity = np.array(UpParity.parity_interest)
        up = np.array(UpParity.up_interest)
        parity_PAT = np.array(UpParity.parity_PAT)
        parity_QPD = np.array(UpParity.parity_QPD)
        up_PAT = np.array(UpParity.up_PAT)
        up_QPD = np.array(UpParity.up_QPD)

        ### plot start
        label_font = 18
        tick_font = 16
        legend_font = 14

        fig, axs = plt.subplots(2, figsize=(6, 8),
                                gridspec_kw={'height_ratios': [3, 4],
                                             'hspace': 0.2}
                                )

        axs[0].errorbar(Up_1D_Time, Up_1D_P1*100, yerr=Up_1D_Std/np.sqrt(20), fmt='o', c='k')
        axs[0].plot(Up_1D_Time, Up_1D_P1Fit*100, c='g')
        axs[0].set_xlim([3.8, 15.2])
        axs[0].set_ylim([1.18, 2.01])
        axs[0].set_yticks([1.2, 1.4, 1.6, 1.8, 2.0])

        axs[0].set_ylabel('P$_{1}$ (%)', fontsize=label_font, family='Arial')
        axs[0].tick_params(labelsize=tick_font)
        f_l = 100
        f_r = 620
        loc = 4
        axs[1].plot(f_interest, parity, 'r', linewidth=3, label='             ')
        axs[1].plot(f_interest, up, 'b', linewidth=3, label='             ')
        axs[1].plot(f_interest, up_PAT, 'b--', linewidth=3, label='             ')
        axs[1].set_ylim([8e1, 1e4])
        axs[1].set_yscale('log')
        axs[1].set_ylabel('Rate ($\mathrm{s}^{-1}$)', fontsize=label_font, family='Arial')
        axs[1].tick_params(labelsize=tick_font)
        axs[1].legend(loc=loc, prop={'size': 12})
        axs[1].set_xlim([f_l, f_r])
        axs[1].set_xlabel('Frequency (GHz)', fontsize=label_font, family='Arial')
        fig.align_ylabels(axs)
        plt.tight_layout()
        path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        plt.savefig(path + '\ParityUpCorrelation.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()

