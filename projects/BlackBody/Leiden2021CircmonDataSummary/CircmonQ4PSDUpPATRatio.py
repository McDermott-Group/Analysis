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
        label_font = 14
        tick_font = 12
        legend_font = 14
        # plt.rcParams['text.usetex'] = True


        plt.figure(figsize=(6, 4))

        f_l = 100
        f_r = 620

        plt.plot(f_interest, np.divide(parity_PAT, parity), 'r', linewidth=3, label='$\Gamma_{\mathrm{P}}$')
        # axs[2].plot(f_interest, np.divide(parity_QPD, parity), 'r--', linewidth=3, label='$\Gamma_{p, QPD}$ ratio')
        plt.plot(f_interest, np.divide(up_PAT, up), 'b', linewidth=3, label='$\Gamma_{\uparrow}$')
        # axs[2].plot(f_interest, np.divide(up_QPD, up), 'b--', linewidth=3, label='$\Gamma_{01, QPD}$ ratio')

        plt.xlim([f_l, f_r])
        plt.ylim([0.1, 1])
        plt.xlabel('Transmitter frequency (GHz)', fontsize=label_font)
        plt.ylabel('PAT1/PAT', fontsize=label_font)
        plt.tick_params(labelsize=tick_font)
        plt.legend(loc=1, prop={'size': 15})
        plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        plt.tick_params(axis="x", direction="in", which='both')
        plt.tick_params(axis="y", direction="in", which='both')
        plt.tick_params(axis="x", width=1, length=6, which='both')
        plt.tick_params(axis="y", width=1, length=6, which='both')

        # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
        # plt.savefig(path + '\S_PATRatio.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()

