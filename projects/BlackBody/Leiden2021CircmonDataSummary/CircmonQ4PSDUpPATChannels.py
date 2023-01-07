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
        # Q4_PSD_Aug = np.loadtxt(Q4_PSD_file_Aug, skiprows=20)

        Q4_P1_file_Aug = "2021AugCircRadiator_Q4_P1_Data.txt"
        # Q4_P1_Aug = np.loadtxt(Q4_P1_file_Aug, skiprows=20)

        Q4_Up_file_Aug = "2021AugCircRadiator_Q4_Up_Data.txt"
        # Q4_Up_Aug = np.loadtxt(Q4_Up_file_Aug, skiprows=20)

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
        label_font = 14
        tick_font = 12
        legend_font = 14
        # plt.rcParams['text.usetex'] = True


        plt.figure(figsize=(6, 4))

        f_l = 100
        f_r = 620

        plt.plot(f_interest, parity_PAT, 'r', linewidth=3, label='                    ')
        plt.plot(f_interest, up_PAT, 'b', linewidth=3, label='                    ')
        plt.plot(f_interest, parity_QPD, 'k', linewidth=3, label='                    ')

        plt.xlim([f_l, f_r])
        plt.xlabel('Transmitter frequency (GHz)', fontsize=label_font)
        plt.yscale('log')
        plt.ylim([1e1, 1e4])
        plt.ylabel('Rate ($\mathrm{s}^{-1}$)', fontsize=label_font, family='Arial')
        # plt.ylabel('PAT1/PAT', fontsize=label_font)
        plt.tick_params(labelsize=tick_font)
        plt.legend(loc=4, prop={'size': 14}, frameon=False)
        # plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        plt.tick_params(axis="x", direction="in", which='both')
        plt.tick_params(axis="y", direction="in", which='both')
        plt.tick_params(axis="x", width=1, length=6, which='both')
        plt.tick_params(axis="y", width=1, length=3, which='minor')
        plt.tick_params(axis="y", width=1, length=6, which='major')



        path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        plt.savefig(path + '\S_PATChannels.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()

