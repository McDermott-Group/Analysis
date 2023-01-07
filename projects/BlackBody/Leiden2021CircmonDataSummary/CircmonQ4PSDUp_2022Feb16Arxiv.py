import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, UpAndParity, getPhotonRate
import numpy as np
import matplotlib.font_manager as font_manager
from scipy import interpolate
from scipy.integrate import quad

if 0:
    """
    Import CST files and Junction parameters
    """
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff_circ = 75 * 1e-21  # Commonly used (50-100)
    C_eff_SFQ = 75 * 1e-21  # Commonly used (50-100)

    JCirc = [4.8 * 1e3, None, 0, 320 * 123 * 2, "Radiator"]  # [R, L, C, A]
    JSFQ_weak = [16 * 1e3, None, 0, 100 * 200, "Radiator"]  # [R, L, C, A]
    JSFQ_strong = [8 * 1e3, None, 0, 100 * 200 * 2, "Radiator"]  # [R, L, C, A]

    JQ1 = [16.6 * 1e3, 18.3 * 1e-9, 0, 193.8 * 121.8, "Receiver"]  #
    JQ2 = [13.2 * 1e3, 14.6 * 1e-9, 0, 350 * 126.6, "Receiver"]  #
    JQ3 = [19 * 1e3, 21 * 1e-9, 0, 310 * 126, "Receiver"]  # some issue with Q3
    JQ4 = [15 * 1e3, 19.9 * 1e-9, 0, 184.4 * 122.5 * 2, "Receiver"]  #

    ### up to 500 GHz
    fileCirc = "2021Aug_J2.txt"
    fileQ1Aug = "2021Aug_Q1_leads.txt"
    fileQ2Aug = "2021Aug_Q2.txt"
    fileQ3Aug = "2021Aug_Q3.txt"
    fileQ4Aug = "2021Aug_Q4.txt"

    ### up to 1000 GHz
    # fileSFQ = "testpad_1.5THz.txt"
    fileSFQ_4 = "SFQ_4rectmons.txt"
    fileSFQ = "SFQ_1THz.txt"
    # fileQ_l = "Q_l_with-leads_1.5THz.txt"
    fileQ1 = "Q1_full-chip.txt"
    fileQ2 = "Q2_full-chip.txt"
    fileQ3 = "Q3_full-chip.txt"
    fileQ4 = "Q4_full-chip.txt"

    Circ = AntennaCoupling()
    Circ.import_data(fileCirc, JCirc, C_eff=C_eff_circ)
    f_Circ = Circ.Antenna["f"]
    ecCirc = Circ.Antenna["e_c_dB"]
    eCirc = Circ.Antenna["e_c"]

    SFQ_weak = AntennaCoupling()
    SFQ_weak.import_data(fileSFQ, JSFQ_weak, C_eff=C_eff_SFQ)
    f_SFQ = SFQ_weak.Antenna["f"]
    ecSFQ_weak = SFQ_weak.Antenna["e_c_dB"]
    eSFQ_weak = SFQ_weak.Antenna["e_c"]

    SFQ_strong = AntennaCoupling()
    SFQ_strong.import_data(fileSFQ, JSFQ_strong, C_eff=C_eff_SFQ)
    ecSFQ_strong = SFQ_strong.Antenna["e_c_dB"]
    eSFQ_strong = SFQ_strong.Antenna["e_c"]

    SFQ_4 = AntennaCoupling()
    SFQ_4.import_data(fileSFQ_4, JSFQ_weak, C_eff=C_eff_SFQ)
    ecSFQ_4 = SFQ_4.Antenna["e_c_dB"]
    eSFQ_4 = SFQ_4.Antenna["e_c"]
    PhotonFlux = SFQ_4.Al_Wall["PhotonFlux"]

    Q1 = AntennaCoupling()
    Q1.import_data(fileQ1, JQ1, C_eff=C_eff_circ)
    eQ1 = Q1.Antenna["e_c"]
    ecQ1 = Q1.Antenna["e_c_dB"]
    AreaQ1 = Q1.Receiver["Area"]

    Q2 = AntennaCoupling()
    Q2.import_data(fileQ2, JQ2, C_eff=C_eff_circ)
    eQ2 = Q2.Antenna["e_c"]
    ecQ2 = Q2.Antenna["e_c_dB"]
    AreaQ2 = Q2.Receiver["Area"]

    Q4 = AntennaCoupling()
    Q4.import_data(fileQ4, JQ4, C_eff=C_eff_circ)
    eQ4 = Q4.Antenna["e_c"]
    ecQ4 = Q4.Antenna["e_c_dB"]
    AreaQ4 = Q4.Receiver["Area"]

    f_scale = np.sqrt(e_eff / 6.0)
    f_Circ = f_Circ * f_scale / 1e9
    f_SFQ = f_SFQ * f_scale / 1e9

    if 0:  # radiator
        plt.plot(f_Circ, eCirc, label='circ')
        plt.plot(f_SFQ, eSFQ_weak, label='weak single')
        plt.plot(f_SFQ, eSFQ_strong, label='strong single')
        plt.plot(f_SFQ, eSFQ_4, '--', label='weak full')
        plt.ylim([3e-3, 5e-1])
    if 0:  # receiver
        plt.plot(f_SFQ, eQ1, c='r', label='Q1')
        plt.plot(f_SFQ, eQ2, c='b', label='Q2')
        plt.plot(f_SFQ, eQ4, c='k', label='Q4')
        plt.ylim([5e-3, 8e-2])

    # plt.yscale('log')
    # plt.xlim([50, 500])
    # # plt.ylim([1e-3, 5e-1])
    # plt.legend()
    # plt.grid(which='both')
    # plt.show()

    Tbb4 = 542.3e-3
    PRQ4 = getPhotonRate(eQ4, f_SFQ, Tbb4)
    PR, f_used = PRQ4[1], PRQ4[2]
    print('PRQ4=', PRQ4[0])

    # UpParity = UpAndParity()
    # freq_data = f_used
    # parity_data = PR
    # UpParity.import_data(freq_data, parity_data)
    # UpRate = UpParity.UpRate
    # print('UpRate=', sum(UpRate))

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
        label_font = 18
        tick_font = 16
        legend_font = 14
        # plt.rcParams['text.usetex'] = True

        legend_font = font_manager.FontProperties(
                                           weight='bold',
                                           style='normal', size=legend_font)

        fig, axs = plt.subplots(2, figsize=(6, 8),
                                gridspec_kw={'height_ratios': [3, 4],
                                             'hspace': 0.2}
                                )

        axs[0].errorbar(Up_1D_Time, Up_1D_P1, yerr=Up_1D_Std/np.sqrt(20), fmt='o', c='k')
        axs[0].plot(Up_1D_Time, Up_1D_P1Fit, c='g')
        axs[0].set_xlim([3.8, 15.2])
        axs[0].set_ylim([0.0118, 0.0205])
        axs[0].set_yticks([0.012, 0.016, 0.020])

        # axs[0].set_xlabel(r"$\Delta t$ ($\mu$s)", fontsize=label_font)
        # axs[0].set_xlabel("\n", fontsize=label_font)
        axs[0].set_ylabel('P($M_{1}(0)|M_{2}(1)$)', fontsize=label_font, family='Arial')
        axs[0].tick_params(labelsize=tick_font)
        # axs[0].legend(loc=loc, fontsize=legend_font)

        f_l = 80
        f_r = 650
        loc = 4
        axs[1].plot(parity_freq_data, parity_data, 'r', linewidth=3, label='$\Gamma_{P, \mathrm{tot}}$')
        # axs[1].plot(f_interest, parity_PAT, 'r--', linewidth=3, label='$\Gamma_{P, \mathrm{PAT1}}$')
        axs[1].plot(up_freq_data, up_data, 'b', linewidth=3, label='$\Gamma_{01, \mathrm{tot}}$')
        # axs[1].plot(f_interest, parity_QPD, 'r:', linewidth=3, label='$\Gamma_{p, QPD}$')
        axs[1].plot(f_interest, up_PAT, 'b--', linewidth=3, label='$\Gamma_{01, \mathrm{PAT1}}$')
        axs[1].set_ylim([8e1, 1e4])
        axs[1].set_yscale('log')
        axs[1].set_ylabel('Rate ($\mathrm{s}^{-1}$)', fontsize=label_font, family='Arial')
        axs[1].tick_params(labelsize=tick_font)
        axs[1].legend(loc=loc, prop=legend_font)
        axs[1].set_xlim([f_l, f_r])
        axs[1].set_xlabel('Frequency (GHz)', fontsize=label_font, family='Arial')

        # axs[2].plot(f_interest, np.divide(parity_PAT, parity), 'r:', linewidth=3, label='$\Gamma_{P, \mathrm{PAT1}}/\Gamma_{P, \mathrm{PAT}}$')
        # # axs[2].plot(f_interest, np.divide(parity_QPD, parity), 'r--', linewidth=3, label='$\Gamma_{p, QPD}$ ratio')
        # axs[2].plot(f_interest, np.divide(up_PAT, up), 'b:', linewidth=3, label='$\Gamma_{01, \mathrm{PAT1}}/\Gamma_{01, \mathrm{PAT}}$')
        # # axs[2].plot(f_interest, np.divide(up_QPD, up), 'b--', linewidth=3, label='$\Gamma_{01, QPD}$ ratio')
        # axs[2].set_xlim([f_l, f_r])
        # axs[2].set_ylim([0, 1])
        # axs[2].set_xlabel('Frequency (GHz)', fontsize=label_font, family='Arial')
        # axs[2].set_ylabel('Primary Ratio', fontsize=label_font, family='Arial')
        # axs[2].tick_params(labelsize=tick_font)
        # axs[2].legend(loc=loc, prop=legend_font)

        fig.align_ylabels(axs)

        plt.tight_layout()
        # plt.title('Up rate baseline hand chosen 200 Hz')
        # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        # plt.savefig(path + '\ParityUpCorrelation.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()

        # fig, axs = plt.subplots(2, sharex='col', figsize=(6, 6),
        #                         gridspec_kw={'height_ratios': [4, 3],
        #                                      'hspace': 0.1}
        #                         )
        #
        # f_l = 80
        # f_r = 650
        # loc = 4
        # axs[0].plot(parity_freq_data, parity_data, 'r', linewidth=3, label='$\Gamma_{p, Meas}$')
        # axs[0].plot(f_interest, parity_PAT, 'r--', linewidth=3, label='$\Gamma_{p, PAT}$')
        # axs[0].plot(up_freq_data, up_data, 'b', linewidth=3, label='$\Gamma_{01, Meas}$')
        #
        #
        # # axs[1].plot(f_interest, parity_QPD, 'r:', linewidth=3, label='$\Gamma_{p, QPD}$')
        # axs[0].plot(f_interest, up_PAT, 'b--', linewidth=3, label='$\Gamma_{01, PAT}$')
        # # axs[1].plot(f_interest, up_QPD, 'b:', linewidth=3, label='$\Gamma_{01, QPD}=\Gamma_{p, QPD}$')
        # # axs[1].set_ylim([8e1, 1e4])
        # # axs[1].set_yscale('log')
        # # axs[1].set_ylabel('Rate (Hz)')
        # # axs[1].legend(loc=loc)
        #
        # axs[0].set_ylim([8e1, 1e4])
        # axs[0].set_yscale('log')
        # axs[0].set_ylabel('Rate (Hz)', fontsize=label_font)
        # axs[0].tick_params(labelsize=tick_font)
        # axs[0].legend(loc=loc, fontsize=legend_font)
        #
        # axs[1].plot(f_interest, np.divide(parity_PAT, parity), 'r:', linewidth=3, label='$\Gamma_{p, PAT}/\Gamma_{p, J}$')
        # # axs[2].plot(f_interest, np.divide(parity_QPD, parity), 'r--', linewidth=3, label='$\Gamma_{p, QPD}$ ratio')
        # axs[1].plot(f_interest, np.divide(up_PAT, up), 'b:', linewidth=3, label='$\Gamma_{01, PAT}/\Gamma_{01, J}$')
        # # axs[2].plot(f_interest, np.divide(up_QPD, up), 'b--', linewidth=3, label='$\Gamma_{01, QPD}$ ratio')
        #
        # axs[1].set_xlim([f_l, f_r])
        # axs[1].set_ylim([0, 1])
        # axs[1].set_xlabel('Frequency (GHz)', fontsize=label_font)
        # axs[1].set_ylabel('Ratio', fontsize=label_font)
        # axs[1].tick_params(labelsize=tick_font)
        # axs[1].legend(loc=loc, fontsize=legend_font)
        #
        # # plt.title('Up rate baseline hand chosen 200 Hz')
        # # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        # # plt.savefig(path + '\ParityUpCorrelation.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        # plt.show()



###Delta = 46 GHz
##S-(245GHz) = 4.2387
##S+(245GHz) = 6.0544
