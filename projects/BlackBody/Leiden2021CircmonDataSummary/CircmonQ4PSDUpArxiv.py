import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, UpAndParity, getPhotonRate
import numpy as np
from scipy import interpolate
from scipy.integrate import quad

if 1:
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
        Q1_PSD_file_Aug = "2021AugCircRadiator_Q1_PSD_Data.txt"
        Q2_PSD_file_Aug = "2021AugCircRadiator_Q2_PSD_Data.txt"
        Q4_PSD_file_Aug = "2021AugCircRadiator_Q4_PSD_Data.txt"
        Q1_PSD_Aug = np.loadtxt(Q1_PSD_file_Aug, skiprows=20)
        Q2_PSD_Aug = np.loadtxt(Q2_PSD_file_Aug, skiprows=20)
        Q4_PSD_Aug = np.loadtxt(Q4_PSD_file_Aug, skiprows=20)

        Q1_P1_file_Aug = "2021AugCircRadiator_Q1_P1_Data.txt"
        Q2_P1_file_Aug = "2021AugCircRadiator_Q2_P1_Data.txt"
        Q3_P1_file_Aug = "2021AugCircRadiator_Q3_P1_Data.txt"
        Q4_P1_file_Aug = "2021AugCircRadiator_Q4_P1_Data.txt"
        Q1_P1_Aug = np.loadtxt(Q1_P1_file_Aug, skiprows=20)
        Q2_P1_Aug = np.loadtxt(Q2_P1_file_Aug, skiprows=20)
        Q3_P1_Aug = np.loadtxt(Q3_P1_file_Aug, skiprows=20)
        Q4_P1_Aug = np.loadtxt(Q4_P1_file_Aug, skiprows=20)

        Q1_Up_file_Aug = "2021AugCircRadiator_Q1_Up_Data.txt"
        Q2_Up_file_Aug = "2021AugCircRadiator_Q2_Up_Data.txt"
        Q4_Up_file_Aug = "2021AugCircRadiator_Q4_Up_Data.txt"
        Q1_Up_Aug = np.loadtxt(Q1_Up_file_Aug, skiprows=20)
        Q2_Up_Aug = np.loadtxt(Q2_Up_file_Aug, skiprows=20)
        Q4_Up_Aug = np.loadtxt(Q4_Up_file_Aug, skiprows=20)

        S_param_file = "S_minusplus.txt"
        # S_param = np.loadtxt(S_param_file)
        # plt.plot(S_param[:, 0], S_param[:, 1])
        # plt.plot(S_param[:, 0], S_param[:, 2])
        # plt.xlim([2, 18])
        # plt.ylim([0, 18])
        # plt.show()

        """Working on the theory"""
        UpParity = UpAndParity()

        # parity_freq_data = Q4_PSD_Aug[:, 0] * f_DAC_Aug
        # parity_data = Q4_PSD_Aug[:, 1]
        # up_freq_data = Q4_Up_Aug[:, 0] * f_DAC_Aug
        # up_data = Q4_Up_Aug[:, 1]

        UpParity.import_data(Q4_PSD_file_Aug, Q4_Up_file_Aug, S_param_file)
        # UpRate = np.array(UpParity.UpRate)
        # X_QP = np.array(UpParity.X_QP)

        BB_UpRate = 52.7  # Hz, calculated from effective BB temp and Yale's theory

        # plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, Q1_PSD_Aug[:, 1])
        # plt.plot(Q2_PSD_Aug[:, 0]*f_DAC_Aug, Q2_PSD_Aug[:, 1])
        # plt.plot(Q4_PSD_Aug[:, 0]*f_DAC_Aug, Q4_PSD_Aug[:, 1], label='Parity')

        # plt.plot(Q1_P1_Aug[:, 0]*f_DAC_Aug, Q1_P1_Aug[:, 1])
        # plt.plot(Q2_P1_Aug[:, 0]*f_DAC_Aug, Q2_P1_Aug[:, 1])
        # plt.plot(Q3_P1_Aug[:, 0]*f_DAC_Aug, Q3_P1_Aug[:, 1])
        # plt.plot(Q4_P1_Aug[:, 0]*f_DAC_Aug, Q4_P1_Aug[:, 1], label='Steady')

        # plt.plot(Q1_Up_Aug[:, 0]*f_DAC_Aug, Q1_Up_Aug[:, 1])
        # plt.plot(Q2_Up_Aug[:, 0]*f_DAC_Aug, Q2_Up_Aug[:, 1])
        # plt.plot(Q4_Up_Aug[:, 0]*f_DAC_Aug, Q4_Up_Aug[:, 1], label='Up')

        ### interpolate since parity rate has higher density
        f = interpolate.interp1d(Q4_PSD_Aug[:, 0], Q4_PSD_Aug[:, 1])
        f_UpParity = interpolate.interp1d(Q4_PSD_Aug[:, 0], UpRate)
        PSD = np.arange(0, 100, 0.1)
        Q4_PSD_Interpolated = f(Q4_Up_Aug[:, 0])
        Q4_Up_Interpolated = f_UpParity(Q4_Up_Aug[:, 0])

        alpha = 8e6

        fig, ax = plt.subplots(1, figsize=(8, 8))
        ax.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, Q4_PSD_Aug[:, 1], 'b', label='$\Gamma_{P}$', linewidth=4)
        # ax[0].plot(PSD * f_DAC_Aug, f(PSD), label='interploted')
        # ax[0].plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, f(Q4_Up_Aug[:, 0]), label='interploted for up rate')
        ax.plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, Q4_Up_Aug[:, 1], 'r', label='$\Gamma_{01}$', linewidth=4)

        # ax.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, alpha * X_QP, 'y', label='X_QP')
        ax.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, UpRate + BB_UpRate, 'g', label='$\Gamma_{01}$ from $\Gamma_{P}$',
                linewidth=4)
        # ax.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, UpRate + alpha * X_QP, 'k', label='Sum')

        # ax[0].set_title('Circmon Q4 Parity vs Up Rate and P1 Steady state')
        # ax_02 = ax.twinx()
        # # ax_02.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, alpha*X_QP, label='X_QP')
        # ax_02.plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, np.divide(Q4_Up_Interpolated, Q4_Up_Aug[:, 1]), 'r--',
        #            label='$\Gamma_{01 from P}/\Gamma_{01}$')
        # ax_02.set_ylabel('$\Gamma_{01 from P}/\Gamma_{01}$', color='red')
        plt.yscale('log')
        plt.xlim([0, 500])
        plt.ylim([1e1, 1e4])

        plt.xlabel('Radiator Freq (GHz)')
        plt.legend()
        plt.show()

        # fig, ax = plt.subplots(2, figsize=(8, 8))
        # ax[0].plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, Q4_PSD_Aug[:, 1], label='$\Gamma_{P}$')
        # ax[0].plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, alpha*X_QP, label='X_QP')
        # # ax[0].plot(PSD * f_DAC_Aug, f(PSD), label='interploted')
        # # ax[0].plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, f(Q4_Up_Aug[:, 0]), label='interploted for up rate')
        # ax[0].plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, Q4_Up_Aug[:, 1], label='$\Gamma_{01}$')
        # ax[0].plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, UpRate, label='$\Gamma_{01}$ from $\Gamma_{P}$')
        # ax[0].set_xlim([0, 550])
        # ax[0].set_ylabel('Rate (Hz)')
        # ax[0].set_yscale('log')
        # ax[0].legend(loc=4)
        # ax[0].set_title('Circmon Q4 Parity vs Up Rate and P1 Steady state')
        #
        # # ax_02 = ax[0].twinx()
        # # ax_02.plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, np.divide(Q4_PSD_Interpolated, Q4_Up_Aug[:, 1]), 'r--',
        # #            label='$\Gamma_{p}/\Gamma_{01}$')
        # # ax_02.set_ylabel('$\Gamma_{p}/\Gamma_{01}$', color='red')
        # # ax_02.set_ylim([1, 7])
        #
        #
        # ax_02 = ax[0].twinx()
        # # ax_02.plot(Q4_PSD_Aug[:, 0] * f_DAC_Aug, alpha*X_QP, label='X_QP')
        # ax_02.plot(Q4_Up_Aug[:, 0] * f_DAC_Aug, np.divide(Q4_Up_Interpolated, Q4_Up_Aug[:, 1]), 'r--',
        #            label='$\Gamma_{01 from P}/\Gamma_{01}$')
        # ax_02.set_ylabel('$\Gamma_{01 from P}/\Gamma_{01}$', color='red')
        # ax_02.set_ylim([0.05, 0.3])

        # ax[1].plot(Q4_P1_Aug[:, 0]*f_DAC_Aug, Q4_P1_Aug[:, 1], label='P1')
        # ax[1].set_ylabel('P1')
        # ax[1].set_yscale('log')
        # ax[1].set_xlim([0, 550])

        # plt.xlim([50, 600])
        # plt.ylim([1e2, 1e4])
        # plt.xlabel('Radiator Freq (GHz)')
        # plt.legend()
        # plt.show()


###Delta = 46 GHz
##S-(245GHz) = 4.2387
##S+(245GHz) = 6.0544
