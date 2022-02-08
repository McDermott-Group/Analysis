import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, getPhotonRate
import numpy as np

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



    JQ1 = [16.6 * 1e3, 18.3 * 1e-9, 0, 193.8 * 121.8 * 2, "Receiver"]  #
    JQ2 = [13.2 * 1e3, 14.6 * 1e-9, 0, 330 * 128.0, "Receiver"]  #
    JQ3 = [19 * 1e3, 21 * 1e-9, 0, 310 * 126, "Receiver"]  # some issue with Q3
    JQ4 = [15 * 1e3, 19.9 * 1e-9, 0, 184.4 * 122.5 * 2, "Receiver"]  #

    ### up to 500 GHz
    fileCirc = "2021Aug_J2.txt"
    fileQ1Aug = "2021Aug_Q1_leads.txt"
    fileQ2Aug = "2021Aug_Q2.txt"
    fileQ3Aug = "2021Aug_Q3.txt"
    fileQ4Aug = "2021Aug_Q4.txt"

    ### up to 1000 GHz
    fileCircTest = "testpad_circmon1THz.txt"
    # fileSFQ = "testpad_1.5THz.txt"
    fileSFQ_4 = "SFQ_4rectmons.txt"
    fileSFQ = "SFQ_1THz.txt"
    # fileQ_l = "Q_l_with-leads_1.5THz.txt"
    fileQ1 = "Q1_full-chip.txt"
    fileQ2 = "Q2_full-chip.txt"
    # fileQ2 = "Q2_overetch.txt"
    # fileQ2 = "Q2_non_overetch.txt"
    fileQ3 = "Q3_full-chip.txt"
    fileQ4 = "Q4_full-chip.txt"

    Circ_test = AntennaCoupling()
    Circ_test.import_data(fileCircTest, JCirc, C_eff=C_eff_circ)
    f_Circ_test = Circ_test.Antenna["f"]
    ecCirc_test = Circ_test.Antenna["e_c_dB"]
    eCirc_test = Circ_test.Antenna["e_c"]
    PhotonFlux_CircTest = Circ_test.Al_Wall["PhotonFlux"]

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
    # Q2.import_data(fileQ2Aug, JQ2, C_eff=C_eff_circ)
    f_Q2 = Q2.Antenna["f"]
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
    f_Circ_test = f_Circ_test * f_scale / 1e9
    f_SFQ = f_SFQ * f_scale / 1e9
    f_Q2 = f_Q2 * f_scale / 1e9

    if 0: # radiator
        plt.plot(f_Circ, eCirc, label='circ')
        plt.plot(f_SFQ, eSFQ_weak, label='weak single')
        plt.plot(f_SFQ, eSFQ_strong, label='strong single')
        plt.plot(f_SFQ, eSFQ_4, '--', label='weak full')
        plt.ylim([3e-3, 5e-1])
    if 0: # receiver
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

        # plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, Q1_PSD_Aug[:, 1])
        # plt.plot(Q2_PSD_Aug[:, 0]*f_DAC_Aug, Q2_PSD_Aug[:, 1])
        # plt.plot(Q4_PSD_Aug[:, 0]*f_DAC_Aug, Q4_PSD_Aug[:, 1])

        # plt.plot(Q1_P1_Aug[:, 0]*f_DAC_Aug, Q1_P1_Aug[:, 1])
        # plt.plot(Q2_P1_Aug[:, 0]*f_DAC_Aug, Q2_P1_Aug[:, 1])
        # plt.plot(Q3_P1_Aug[:, 0]*f_DAC_Aug, Q3_P1_Aug[:, 1])
        # plt.plot(Q4_P1_Aug[:, 0]*f_DAC_Aug, Q4_P1_Aug[:, 1])

        # plt.plot(Q1_Up_Aug[:, 0]*f_DAC_Aug, Q1_Up_Aug[:, 1])
        # plt.plot(Q2_Up_Aug[:, 0]*f_DAC_Aug, Q2_Up_Aug[:, 1])
        # plt.plot(Q4_Up_Aug[:, 0]*f_DAC_Aug, Q4_Up_Aug[:, 1])

    if 1:
        """SepSFQRectmonStrong Data"""
        f_DAC_Sep = 4.604
        Q1_PSD_file_Sep = "2021SepSFQStrongRadiator_Q1_PSD_Data.txt"
        Q2_PSD_file_Sep = "2021SepSFQStrongRadiator_Q2_PSD_Data.txt"
        Q4_PSD_file_Sep = "2021SepSFQStrongRadiator_Q4_PSD_Data.txt"
        Q1_PSD_Sep = np.loadtxt(Q1_PSD_file_Sep, skiprows=20)
        Q2_PSD_Sep = np.loadtxt(Q2_PSD_file_Sep, skiprows=20)
        Q4_PSD_Sep = np.loadtxt(Q4_PSD_file_Sep, skiprows=20)

        # plt.plot(Q1_PSD_Sep[:, 0]*f_DAC_Sep, Q1_PSD_Sep[:, 1])
        # plt.plot(Q2_PSD_Sep[:, 0]*f_DAC_Sep, Q2_PSD_Sep[:, 1])
        # plt.plot(Q4_PSD_Sep[:, 0]*f_DAC_Sep, Q4_PSD_Sep[:, 1])
        # plt.show()

    if 1:
        """OctSFQRectmonWeak Data"""
        f_SIM = 0.9676
        Q1_PSD_file_OctWeak = "2021OctSFQWeakRadiator_Q1_PSD_Data.txt"
        Q2_PSD_file_OctWeak = "2021OctSFQWeakRadiator_Q2_PSD_Data.txt"
        Q4_PSD_file_OctWeak = "2021OctSFQWeakRadiator_Q4_PSD_Data.txt"
        Q1_PSD_OctWeak = np.loadtxt(Q1_PSD_file_OctWeak, skiprows=20)
        Q2_PSD_OctWeak = np.loadtxt(Q2_PSD_file_OctWeak, skiprows=20)
        Q4_PSD_OctWeak = np.loadtxt(Q4_PSD_file_OctWeak, skiprows=20)

        Q1_PSD_file_OctWeak_Wiggles = "2021OctSFQWeakRadiatorWiggles_Q1_PSD_Data.txt"
        Q2_PSD_file_OctWeak_Wiggles = "2021OctSFQWeakRadiatorWiggles_Q2_PSD_Data.txt"
        Q4_PSD_file_OctWeak_Wiggles = "2021OctSFQWeakRadiatorWiggles_Q4_PSD_Data.txt"
        Q1_PSD_OctWeak_Wiggles = np.loadtxt(Q1_PSD_file_OctWeak_Wiggles, skiprows=20)
        Q2_PSD_OctWeak_Wiggles = np.loadtxt(Q2_PSD_file_OctWeak_Wiggles, skiprows=20)
        Q4_PSD_OctWeak_Wiggles = np.loadtxt(Q4_PSD_file_OctWeak_Wiggles, skiprows=20)

        Q1_P1_file_OctWeak = "2021OctSFQWeakRadiator_Q1_P1_Data.txt"
        Q2_P1_file_OctWeak = "2021OctSFQWeakRadiator_Q2_P1_Data.txt"
        Q4_P1_file_OctWeak = "2021OctSFQWeakRadiator_Q4_P1_Data.txt"
        Q1_P1_OctWeak = np.loadtxt(Q1_P1_file_OctWeak, skiprows=20)
        Q2_P1_OctWeak = np.loadtxt(Q2_P1_file_OctWeak, skiprows=20)
        Q4_P1_OctWeak = np.loadtxt(Q4_P1_file_OctWeak, skiprows=20)

        # plt.plot(Q1_PSD_OctWeak[:, 0]*f_SIM, Q1_PSD_OctWeak[:, 1])
        # plt.plot(Q2_PSD_OctWeak[:, 0]*f_SIM, Q2_PSD_OctWeak[:, 1])
        # plt.plot(Q4_PSD_OctWeak[:, 0]*f_SIM, Q4_PSD_OctWeak[:, 1])
        #
        # plt.plot(Q1_PSD_OctWeak_Wiggles[:, 0]*f_SIM, Q1_PSD_OctWeak_Wiggles[:, 1])
        # plt.plot(Q2_PSD_OctWeak_Wiggles[:, 0]*f_SIM, Q2_PSD_OctWeak_Wiggles[:, 1])
        # plt.plot(Q4_PSD_OctWeak_Wiggles[:, 0]*f_SIM, Q4_PSD_OctWeak_Wiggles[:, 1])

        # plt.plot(Q1_P1_OctWeak[:, 0]*f_SIM, Q1_P1_OctWeak[:, 1])
        # plt.plot(Q2_P1_OctWeak[:, 0]*f_SIM, Q2_P1_OctWeak[:, 1])
        # plt.plot(Q4_P1_OctWeak[:, 0]*f_SIM, Q4_P1_OctWeak[:, 1])

    if 1:
        """OctSFQRectmonStrong Data"""
        f_DAC_Oct = 4.604
        f_SIM = 0.9676
        Q1_PSD_file_OctStrong = "2021OctSFQStrongRadiator_Q1_PSD_Data.txt"
        Q2_PSD_file_OctStrong = "2021OctSFQStrongRadiator_Q2_PSD_Data.txt"
        Q4_PSD_file_OctStrong = "2021OctSFQStrongRadiator_Q4_PSD_Data.txt"
        Q1_PSD_OctStrong = np.loadtxt(Q1_PSD_file_OctStrong, skiprows=20)
        Q2_PSD_OctStrong = np.loadtxt(Q2_PSD_file_OctStrong, skiprows=20)
        Q4_PSD_OctStrong = np.loadtxt(Q4_PSD_file_OctStrong, skiprows=20)

        Q1_P1_file_OctStrong = "2021OctSFQStrongRadiator_Q1_P1_Data.txt"
        Q2_P1_file_OctStrong = "2021OctSFQStrongRadiator_Q2_P1_Data.txt"
        Q4_P1_file_OctStrong = "2021OctSFQStrongRadiator_Q4_P1_Data.txt"
        Q1_P1_OctStrong = np.loadtxt(Q1_P1_file_OctStrong, skiprows=20)
        Q2_P1_OctStrong = np.loadtxt(Q2_P1_file_OctStrong, skiprows=20)
        Q4_P1_OctStrong = np.loadtxt(Q4_P1_file_OctStrong, skiprows=20)

        # plt.plot(Q1_PSD_OctStrong[:, 0]*f_DAC_Oct, Q1_PSD_OctStrong[:, 1])
        # plt.plot(Q2_PSD_OctStrong[:, 0]*f_DAC_Oct, Q2_PSD_OctStrong[:, 1])
        # plt.plot(Q4_PSD_OctStrong[:, 0]*f_DAC_Oct, Q4_PSD_OctStrong[:, 1])

        # plt.plot(Q1_P1_OctStrong[:, 0]*f_SIM, Q1_P1_OctStrong[:, 1])
        # plt.plot(Q2_P1_OctStrong[:, 0]*f_SIM, Q2_P1_OctStrong[:, 1])
        # plt.plot(Q4_P1_OctStrong[:, 0]*f_SIM, Q4_P1_OctStrong[:, 1])

    if 0: # plot P1 data
        f_SIM_OctW = 0.96758 * 1.025 # Xmon match data
        SIM_OctW_Offset = 2 # Xmon match data   2mV is OK or max
        f_SIM_Octs = 0.96758 * 1 # Xmon match data
        plt.plot(Q1_P1_Aug[:, 0]*f_DAC_Aug, Q1_P1_Aug[:, 1], 'r--', label='Q1 P1 Aug')
        plt.plot(Q2_P1_Aug[:, 0]*f_DAC_Aug, Q2_P1_Aug[:, 1], 'b--', label='Q2 P1 Aug')
        plt.plot(Q4_P1_Aug[:, 0]*f_DAC_Aug, Q4_P1_Aug[:, 1], 'k--', label='Q4 P1 Aug')

        plt.plot((Q1_P1_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, 2*Q1_P1_OctWeak[:, 1], 'r-', label='Q1 P1 Oct Weak')
        # plt.plot(Q2_P1_OctWeak[:, 0]*f_SIM_OctW, Q2_P1_OctWeak[:, 1], 'b-', label='Q2 P1 Oct Weak')
        plt.plot((Q4_P1_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, 2*Q4_P1_OctWeak[:, 1], 'k-', label='Q4 P1 Oct Weak')

        plt.plot(Q1_P1_OctStrong[:, 0]*f_SIM_Octs, Q1_P1_OctStrong[:, 1], 'r-.', label='Q1 P1 Oct Strong')
        # plt.plot(Q2_P1_OctStrong[:, 0]*f_SIM_Octs, Q2_P1_OctStrong[:, 1], 'b-.', label='Q1 P1 Oct Strong')
        plt.plot(Q4_P1_OctStrong[:, 0]*f_SIM_Octs, Q4_P1_OctStrong[:, 1], 'k-.', label='Q1 P1 Oct Strong')

    if 0:   # plot PSD data
        f_SIM_OctW = 0.96758 * 1.025 # Xmon match data
        SIM_OctW_Offset = 2 # Xmon match data   2mV is OK or max
        # f_SIM_Octs = 0.96758 * 1 # Xmon match data
        f_DAC_Oct = 4.6   # 4.604
        DAC_Oct_Offset = 0.5    # 0.9
        plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, Q1_PSD_Aug[:, 1], 'r--', label='Q1 PSD Aug')
        plt.plot(Q2_PSD_Aug[:, 0]*f_DAC_Aug, Q2_PSD_Aug[:, 1], 'b--', label='Q2 PSD Aug')
        plt.plot(Q4_PSD_Aug[:, 0]*f_DAC_Aug, Q4_PSD_Aug[:, 1], 'k--', label='Q4 PSD Aug')

        plt.plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', label='Q1 PSD Oct Weak')
        plt.plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1], 'b-', label='Q2 PSD Oct Weak')
        plt.plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1], 'k-', label='Q4 PSD Oct Weak')

        plt.plot((Q1_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q1_PSD_OctStrong[:, 1], 'r-.', label='Q1 PSD Oct Strong')
        plt.plot((Q2_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q2_PSD_OctStrong[:, 1], 'b-.', label='Q2 PSD Oct Strong')
        plt.plot((Q4_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q4_PSD_OctStrong[:, 1], 'k-.', label='Q4 PSD Oct Strong')

        plt.xlim([50, 600])
        # plt.xlim([1e2, 1e5])
        plt.yscale('log')
        plt.legend()
        plt.show()

    if 0:   # parity data cross talk for Aug
        f_SIM_OctW = 0.96758 * 1.025 # Xmon match data
        SIM_OctW_Offset = 2 # Xmon match data   2mV is OK or max
        # f_SIM_Octs = 0.96758 * 1 # Xmon match data
        f_DAC_Oct = 4.64   # 4.604
        DAC_Oct_Offset = 0.5    # 0.9

        k41 = 0     # [, 0.3 (baseline)], [,~0.02(-baseline)]
        k21 = 0     # [, 0.025 (baseline)], [,~0.003(-baseline)]
        k24 = 0     #

        Q1_base = np.mean(Q1_PSD_Aug[:, 1][:8])
        Q2_base = np.mean(Q2_PSD_Aug[:, 1][:8])
        Q4_base = np.mean(Q4_PSD_Aug[:, 1][:8])

        # plt.axhline(y=Q1_base)
        # plt.axhline(y=Q2_base)
        # plt.axhline(y=Q4_base)

        plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, Q1_PSD_Aug[:, 1], 'r--', label='Q1 PSD Aug')
        plt.plot(Q2_PSD_Aug[:, 0]*f_DAC_Aug, Q2_PSD_Aug[:, 1], 'b--', label='Q2 PSD Aug')
        plt.plot(Q4_PSD_Aug[:, 0]*f_DAC_Aug, Q4_PSD_Aug[:, 1], 'k--', label='Q4 PSD Aug')

        # plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, Q1_PSD_Aug[:, 1]-Q1_base, 'r--', label='Q1 PSD Aug')
        # plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, 0.002*(Q1_PSD_Aug[:, 1]-Q1_base), 'r:', label='Q1 ')
        # plt.plot(Q1_PSD_Aug[:, 0]*f_DAC_Aug, 0.005*(Q1_PSD_Aug[:, 1]-Q1_base), 'r:', label='Q1 ')
        # plt.plot(Q2_PSD_Aug[:, 0]*f_DAC_Aug, Q2_PSD_Aug[:, 1]-Q2_base, 'b--', label='Q2 PSD Aug')
        # plt.plot(Q4_PSD_Aug[:, 0]*f_DAC_Aug, Q4_PSD_Aug[:, 1]-Q4_base, 'k--', label='Q4 PSD Aug')
        plt.ylim([1e0, 1e5])

    if 0:   # parity data cross talk for Sep
        f_SIM_OctW = 0.96758 * 1.025 # Xmon match data
        SIM_OctW_Offset = 2 # Xmon match data   2mV is OK or max
        # f_SIM_Octs = 0.96758 * 1 # Xmon match data
        f_DAC_Oct = 4.64   # 4.604
        DAC_Oct_Offset = 0.5    # 0.9

        k41 = 0     # [, 0.33 (baseline)], [, 0.04(-baseline)]
        k21 = 0     # [, 0.025 (baseline)], [, 0.015(-baseline)]
        k24 = 0     #

        Q1_base = np.mean(Q1_PSD_Sep[:, 1][:8])
        Q2_base = np.mean(Q2_PSD_Sep[:, 1][:8])
        Q4_base = np.mean(Q4_PSD_Sep[:, 1][:8])

        # plt.axhline(y=Q1_base)
        # plt.axhline(y=Q2_base)
        # plt.axhline(y=Q4_base)

        # plt.plot(Q1_PSD_Sep[:, 0]*f_DAC_Sep, Q1_PSD_Sep[:, 1], 'r-', label='Q1 PSD Sep')
        # plt.plot(Q2_PSD_Sep[:, 0]*f_DAC_Sep, Q2_PSD_Sep[:, 1], 'b-', label='Q1 PSD Sep')
        # plt.plot(Q4_PSD_Sep[:, 0]*f_DAC_Sep, Q4_PSD_Sep[:, 1], 'k-', label='Q1 PSD Sep')

        # plt.plot(Q1_PSD_Sep[:, 0]*f_DAC_Sep, (Q1_PSD_Sep[:, 1]-Q1_base), 'r-', label='Q1 PSD Sep')
        plt.plot(Q1_PSD_Sep[:, 0]*f_DAC_Sep, 0.014*(Q1_PSD_Sep[:, 1]-Q1_base), 'r-', label='Q1 PSD Sep')
        plt.plot(Q1_PSD_Sep[:, 0]*f_DAC_Sep, 0.016*(Q1_PSD_Sep[:, 1]-Q1_base), 'r-', label='Q1 PSD Sep')
        plt.plot(Q2_PSD_Sep[:, 0]*f_DAC_Sep, Q2_PSD_Sep[:, 1]-Q2_base, 'b-', label='Q1 PSD Sep')
        # plt.plot(Q4_PSD_Sep[:, 0]*f_DAC_Sep, Q4_PSD_Sep[:, 1]-Q4_base, 'k-', label='Q1 PSD Sep')

    if 0:  # parity data cross talk for Oct weak
        f_SIM_OctW = 0.96758 * 1.025  # Xmon match data
        SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
        # f_SIM_Octs = 0.96758 * 1 # Xmon match data
        f_DAC_Oct = 4.64  # 4.604
        DAC_Oct_Offset = 0.5  # 0.9

        k41 = 0     # [, 0.18 (baseline)], [0, 0.12 (-baseline)]
        k21 = 0     # [, 0.012 (baseline)], [0, 0.016 (-baseline)]
        k24 = 0     #
        # Q1_base = 1046.0
        # Q2_base = 12.5
        # Q4_base = 192.0

        Q1_base = np.mean(Q1_PSD_OctWeak[:, 1][:8])
        Q2_base = np.mean(Q2_PSD_OctWeak[:, 1][:8])
        Q4_base = np.mean(Q4_PSD_OctWeak[:, 1][:8])

        plt.axhline(y=Q1_base)
        plt.axhline(y=Q2_base)
        plt.axhline(y=Q4_base)

        plt.plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', label='Q1 PSD Oct Weak')
        plt.plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1], 'b-', label='Q2 PSD Oct Weak')
        plt.plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1], 'k-', label='Q4 PSD Oct Weak')

        # plt.plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, (Q1_PSD_OctWeak[:, 1]-Q1_base), 'r-', label='Q1 PSD Oct Weak')
        # plt.plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, 0.014*(Q1_PSD_OctWeak[:, 1]-Q1_base), 'r:', label='Q1 PSD Oct Weak')
        # plt.plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, 0.016*(Q1_PSD_OctWeak[:, 1]-Q1_base), 'r:', label='Q1 PSD Oct Weak')
        # plt.plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1]-Q2_base, 'b-', label='Q2 PSD Oct Weak')
        # plt.plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1]-Q4_base, 'k-', label='Q4 PSD Oct Weak')

        plt.ylim([1e0, 1e4])

    if 0:  # parity data cross talk for Oct strong
        f_SIM_OctW = 0.96758 * 1.025  # Xmon match data
        SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
        # f_SIM_Octs = 0.96758 * 1 # Xmon match data
        f_DAC_Oct = 4.64  # 4.604
        DAC_Oct_Offset = 0.5  # 0.9

        k41 = 0  # [, 0.18 (baseline)], [ ~(-baseline)]
        k21 = 0  # [, 0.0125 (baseline)], [, 0.01 (-baseline)]
        k24 = 0  #
        # Q1_base = 1016.0
        # Q2_base = 12.8
        # Q4_base = 187.0

        Q1_base = np.mean(Q1_PSD_OctStrong[:, 1][:8])
        Q2_base = np.mean(Q2_PSD_OctStrong[:, 1][:8])
        Q4_base = np.mean(Q4_PSD_OctStrong[:, 1][:8])

        # plt.axhline(y=Q1_base)
        # plt.axhline(y=Q2_base)
        # plt.axhline(y=Q4_base)

        # plt.plot((Q1_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q1_PSD_OctStrong[:, 1], 'r-.', label='Q1 PSD Oct Strong')
        # plt.plot((Q2_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q2_PSD_OctStrong[:, 1], 'b-.', label='Q2 PSD Oct Strong')
        # plt.plot((Q4_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q4_PSD_OctStrong[:, 1], 'k-.', label='Q4 PSD Oct Strong')

        # plt.plot((Q1_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, (Q1_PSD_OctStrong[:, 1]-Q1_base), 'r-.', label='Q1 PSD Oct Strong')
        plt.plot((Q1_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, 0.008*(Q1_PSD_OctStrong[:, 1]-Q1_base), 'r:', label='Q1 PSD Oct Strong')
        plt.plot((Q1_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, 0.01*(Q1_PSD_OctStrong[:, 1]-Q1_base), 'r:', label='Q1 PSD Oct Strong')
        plt.plot((Q2_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q2_PSD_OctStrong[:, 1]-Q2_base, 'b-.', label='Q2 PSD Oct Strong')
        # plt.plot((Q4_PSD_OctStrong[:, 0]+DAC_Oct_Offset)*f_DAC_Oct, Q4_PSD_OctStrong[:, 1]-Q4_base, 'k-.', label='Q4 PSD Oct Strong')

    """
    coefficients
    1. from baseline match, Q1 to Q2 Q4
        Max: k41=0.18, k21=0.012
    2. Data - baseline match, Q1 to Q2 Q4
        Max: k41=~, k21=0.01
    """

    # plt.xlim([50, 600])
    # # plt.xlim([1e2, 1e5])
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

if 0:  # parity data cross talk for Aug
    # f_SIM_OctW = 0.96758 * 1.025  # Xmon match data
    f_SIM_OctW = 0.96758 * 0.95  # Xmon match data
    SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
    f_DAC_Aug = 4.604  # 4.604
    DAC_Oct_Offset = 0.5  # 0.9


    k41 = 0     # [, 0.18 (baseline)], [0, 0.12 (-baseline)]
    k21 = 0     # [, 0.012 (baseline)], [0, 0.016 (-baseline)]
    k24 = 0     #
    # Q1_base = 1046.0
    # Q2_base = 12.5
    # Q4_base = 192.0

    Q1_base = np.mean(Q1_PSD_Aug[:, 1][:8])
    Q2_base = np.mean(Q2_PSD_Aug[:, 1][:8])
    Q4_base = np.mean(Q4_PSD_Aug[:, 1][:8])

    Q1_photon_absorbed = eQ1 * AreaQ1 * PhotonFlux_CircTest/2
    Q2_photon_absorbed = eQ2 * AreaQ2 * PhotonFlux_CircTest/2
    Q4_photon_absorbed = eQ4 * AreaQ4 * PhotonFlux_CircTest/2

    fig, axs = plt.subplots(2, figsize=(8, 8))

    f_l = 50
    f_r = 550

    axs_02 = axs[0].twinx()
    axs_02.plot(f_Circ_test, PhotonFlux_CircTest, color="black", linestyle='--',
                linewidth=4, label='Photon Flux')
    axs_02.set_ylabel("Photon flux $=S/hf$ $(m^{-2} sec^{-1})$", color="black", fontsize=10)
    axs_02.set_xlim([f_l, f_r])
    axs_02.set_ylim([5e12, 5e15])
    axs_02.set_yscale('log')
    axs_02.legend(loc=4)

    axs[0].plot(f_Circ_test, eQ1*AreaQ1, color="red", label='Q1')
    axs[0].plot(f_Circ_test, eQ2*AreaQ2, color="blue", label='Q2')
    axs[0].plot(f_Circ_test, eQ4*AreaQ4, color="black", label='Q4')
    axs[0].set_ylabel('$e_{c}^{r}*A_{eff} (sec^{-1})$', color="black", fontsize=10)
    # axs[0].set_ylabel('$e_{c}$', color="black", fontsize=10)
    axs[0].set_xlim([f_l, f_r])
    axs[0].set_ylim([1e-11, 1e-8])
    axs[0].set_yscale('log')
    axs[0].grid(which='both')
    axs[0].legend(loc=1)

    # axs[0].plot(f_Circ_test, Q1_photon_absorbed, color="red", label='Q1')
    # axs[0].plot(f_Circ_test, Q2_photon_absorbed, color="blue", label='Q2')
    # axs[0].plot(f_Circ_test, Q4_photon_absorbed, color="black", label='Q4')
    # axs[0].set_ylabel('$PAT$', color="black", fontsize=10)
    # axs[0].set_xlim([f_l, f_r])
    # axs[0].set_ylim([3e2, 3e6])
    # axs[0].set_yscale('log')
    # axs[0].grid(which='both')
    # axs[0].legend(loc=1)

    axs[1].plot((Q1_PSD_Aug[:, 0])*f_DAC_Aug, Q1_PSD_Aug[:, 1], 'r-', label='Q1 PSD Aug')
    axs[1].plot((Q2_PSD_Aug[:, 0])*f_DAC_Aug, Q2_PSD_Aug[:, 1], 'b-', label='Q2 PSD Aug')
    axs[1].plot((Q4_PSD_Aug[:, 0])*f_DAC_Aug, Q4_PSD_Aug[:, 1], 'k-', label='Q4 PSD Aug')

    axs[1].set_xlim([f_l, f_r])
    # axs[1].set_ylim([8e2, 8e3]) # Q1
    # axs[1].set_ylim([1e1, 1e3]) # Q2
    # axs[1].set_ylim([1e2, 5e3]) # Q4
    # axs[1].set_ylim([1e1, 6e3])
    axs[1].set_yscale('log')
    axs[1].grid(which='both')
    axs[1].legend(loc=1)
    # plt.yscale('log')
    plt.show()

if 0:  # parity data cross talk for Oct weak
    # f_SIM_OctW = 0.96758 * 1.025  # Xmon match data
    f_SIM_OctW = 0.96758 * 0.95  # Xmon match data
    SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
    f_DAC_Oct = 4.604  # 4.604
    DAC_Oct_Offset = 0.5  # 0.9

    k41 = 0     # [, 0.18 (baseline)], [0, 0.12 (-baseline)]
    k21 = 0     # [, 0.012 (baseline)], [0, 0.016 (-baseline)]
    k24 = 0     #
    # Q1_base = 1046.0
    # Q2_base = 12.5
    # Q4_base = 192.0

    Q1_base = np.mean(Q1_PSD_OctWeak[:, 1][:8])
    Q2_base = np.mean(Q2_PSD_OctWeak[:, 1][:8])
    Q4_base = np.mean(Q4_PSD_OctWeak[:, 1][:8])

    fig, axs = plt.subplots(2, figsize=(8, 8))

    f_l = 50
    f_r = 550

    axs_02 = axs[0].twinx()
    axs_02.plot(f_SFQ, PhotonFlux, color="green", linestyle='--',
                linewidth=4, label='Photon Flux')
    axs_02.set_ylabel("Photon flux $=S/hf$ $(m^{-2} sec^{-1})$", color="black", fontsize=10)
    axs_02.set_xlim([f_l, f_r])
    axs_02.set_ylim([5e12, 5e14])
    axs_02.set_yscale('log')
    axs_02.legend(loc=1)

    axs[0].plot(f_SFQ, eQ1*AreaQ1, color="red", label='Q1')
    axs[0].plot(f_SFQ, eQ2*AreaQ2, color="blue", label='Q2')
    axs[0].plot(f_SFQ, eQ4*AreaQ4, color="black", label='Q4')
    # axs[0].plot(f_SFQ, eSFQ_4, color="purple", label='SFQ Full')
    # axs[0].plot(f_SFQ, eSFQ_weak, color="green", label='SFQ')
    # axs[0].plot(f_SFQ, eQ1*eSFQ_4, color="cyan", label='Tot')
    # axs[0].plot(f_SFQ, eQ2*eSFQ_4, color="cyan", label='Tot')
    # axs[0].plot(f_SFQ, eQ4*eSFQ_4, color="cyan", label='Tot')
    axs[0].set_ylabel('$e_{c}^{r}*A_{eff} (sec^{-1})$', color="black", fontsize=10)
    # axs[0].set_ylabel('$e_{c}$', color="black", fontsize=10)
    axs[0].set_xlim([f_l, f_r])
    axs[0].set_ylim([5e-11, 5e-9])
    # axs[0].set_ylim([5e-3, 6e-2])
    # axs[0].set_ylim([5e-5, 6e-2])   # Q1 tot
    # axs[0].set_ylim([5e-7, 6e-2])   # Q2 tot
    # axs[0].set_ylim([5e-5, 6e-2])   # Q4 tot
    axs[0].set_yscale('log')
    # axs[0].grid(which='both')
    axs[0].legend(loc=2)

    # plt.axhline(y=Q1_base)
    # plt.axhline(y=Q2_base)
    # plt.axhline(y=Q4_base)

    axs[1].plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', label='Q1 PSD Oct Weak')
    axs[1].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1], 'b-', label='Q2 PSD Oct Weak')
    axs[1].plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1], 'k-', label='Q4 PSD Oct Weak')

    axs[1].set_xlim([f_l, f_r])
    # axs[1].set_ylim([8e2, 8e3]) # Q1
    # axs[1].set_ylim([1e1, 1e3]) # Q2
    # axs[1].set_ylim([1e2, 5e3]) # Q4
    axs[1].set_ylim([1e1, 6e3])
    axs[1].set_yscale('log')
    # axs[1].grid(which='both')
    axs[1].legend(loc=1)
    # plt.yscale('log')
    plt.show()

if 0:  # parity data cross talk for Oct weak
    f_SIM_OctW = 0.96758 * 1.025  # Xmon match data
    SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
    f_DAC_Oct = 4.64  # 4.604
    DAC_Oct_Offset = 0.5  # 0.9

    k41 = 0.05     # [, 0.18 (baseline)], [0, 0.12 (-baseline)]
    k21 = 0     # [, 0.012 (baseline)], [0, 0.016 (-baseline)]
    k24 = 0     #
    # Q1_base = 1046.0
    # Q2_base = 12.5
    # Q4_base = 192.0

    Q1_base = np.mean(Q1_PSD_OctWeak[:, 1][:8])
    Q2_base = np.mean(Q2_PSD_OctWeak[:, 1][:8])
    Q4_base = np.mean(Q4_PSD_OctWeak[:, 1][:8])

    Q1_added = Q1_PSD_OctWeak[:, 1]-Q1_base
    Q2_added = Q2_PSD_OctWeak[:, 1]-Q2_base
    Q4_added = Q4_PSD_OctWeak[:, 1]-Q4_base

    Q1_PAT = Q1_added
    Q4_PAT = Q4_added - k41*Q1_PAT
    Q2_PAT = Q2_added - k21*Q1_PAT - k24*Q4_PAT

    Q1_photon_absorbed = eQ1 * AreaQ1 * PhotonFlux/2
    Q2_photon_absorbed = eQ2 * AreaQ2 * PhotonFlux/2
    Q4_photon_absorbed = eQ4 * AreaQ4 * PhotonFlux/2
    etaQ1 = 0.1
    etaQ2 = 0.0
    etaQ4 = 0.4

    fig, axs = plt.subplots(3, figsize=(8, 8))

    f_l = 50
    f_r = 550

    axs_02 = axs[0].twinx()
    axs_02.plot(f_SFQ, PhotonFlux, color="black", linestyle='--',
                linewidth=4, label='Photon Flux')
    axs_02.set_ylabel("Photon flux $=S/hf$ $(m^{-2} sec^{-1})$", color="purple", fontsize=10)
    axs_02.set_xlim([f_l, f_r])
    axs_02.set_ylim([5e12, 5e14])
    axs_02.set_yscale('log')
    axs_02.legend(loc=4)

    # axs[0].plot(f_SFQ, eQ1*AreaQ1, color="red", label='Q1')
    axs[0].plot(f_SFQ, eQ2*AreaQ2, color="blue", label='Q2')
    # axs[0].plot(f_SFQ, eQ4*AreaQ4, color="black", label='Q4')
    # axs[0].plot(f_SFQ, eSFQ_4, color="purple", label='SFQ Full')
    axs[0].set_ylabel('$e_{c}^{r}*A_{eff} (sec^{-1})$', color="black", fontsize=10)
    axs[0].set_xlim([f_l, f_r])
    # axs[0].set_ylim([5e-11, 5e-9])  # Q1
    # axs[0].set_ylim([2e-11, 2e-9])  # Q4
    axs[0].set_ylim([8e-12, 8e-10])  # Q2
    # axs[0].set_ylim([5e-3, 6e-2])
    axs[0].set_yscale('log')
    axs[0].grid(which='both')
    axs[0].legend(loc=1)

    # axs[1].plot(f_SFQ, Q1_photon_absorbed, color="black", label='Q1')
    # axs[1].plot(f_SFQ, Q4_photon_absorbed, color="black", label='Q4')
    axs[1].plot(f_SFQ, Q2_photon_absorbed, color="black", label='Q2')
    axs[1].set_xlim([f_l, f_r])
    # axs[1].set_ylim([1e2, 6e4]) # Q1
    # axs[1].set_ylim([2e2, 2e4]) # Q4
    axs[1].set_ylim([1e1, 1e4]) # Q2
    axs[1].set_yscale('log')
    axs[1].set_ylabel('$e_{c}^{r}*A_{eff}*PhotonFlux/2$', color="black", fontsize=10)
    axs[1].grid(which='both')
    axs[1].legend(loc=1)

    etaQ2 = 0.1
    k21 = 0.01
    k24 = 0.01
    # axs[2].plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1]-Q1_base, 'r-', label='Q1-Base')
    # axs[2].plot(f_SFQ, etaQ1*Q1_photon_absorbed, color="black", label='Q1 Theo')
    # axs[2].plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', label='Q1 Tot')
    # axs[2].plot(f_SFQ, etaQ1*Q1_photon_absorbed + Q1_base, color="black", label='Q1 Theo + Base')
    # axs[2].plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_added, 'b-', label='Q4_added')
    # axs[2].plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_added-k41*Q1_PAT, 'r-', label='Q4_added-Q1_PAT')
    # axs[2].plot(f_SFQ, etaQ4 * Q4_photon_absorbed, color="black", label='Q4 PAT')

    axs[2].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_added, 'k-', label='Q2_added')
    axs[2].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_added-Q1_PAT*k21-Q4_PAT*k24,
                'k-', label='Q2_added-Q1Q4 contribution')
    # axs[2].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PAT*k21, 'r-', label='Q2 from Q1')
    # axs[2].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PAT*k24, 'b-', label='Q2 from Q4')
    axs[2].plot(f_SFQ, etaQ2 * Q2_photon_absorbed, color="purple", label='Q2 PAT')

    axs[2].set_xlim([f_l, f_r])
    # axs[2].set_ylim([1e2, 8e3]) # Q1
    # axs[2].set_ylim([5e1, 8e3]) # Q4
    axs[2].set_ylim([1e0, 6e2]) # Q2
    # axs[1].set_ylim([1e1, 6e3])
    axs[2].set_yscale('log')
    axs[2].grid(which='both')
    axs[2].legend(loc=1)
    # plt.yscale('log')
    plt.show()


### For paper plot

"""
Calculate effective blackbody temperature
"""
# Circmon radiator data
# Tbb1 = 429e-3
# Tbb2 = 522e-3
# Tbb4 = 543e-3
# PRQ1 = getPhotonRate(eQ1, f_SFQ, Tbb1)
# PRQ2 = getPhotonRate(eQ2, f_SFQ, Tbb2)
# PRQ4 = getPhotonRate(eQ4, f_SFQ, Tbb4)
# print('PRQ1=', PRQ1[0])
# print('PRQ2=', PRQ2[0])
# print('PRQ4=', PRQ4[0])

### SFQ weak radiator data
# Tbb1 = 407e-3   # 407,408
# Tbb2 = 459e-3   # 458,460
# Tbb4 = 488e-3   # 487,489
# PRQ1 = getPhotonRate(eQ1, f_SFQ, Tbb1)
# PRQ2 = getPhotonRate(eQ2, f_SFQ, Tbb2)
# PRQ4 = getPhotonRate(eQ4, f_SFQ, Tbb4)
# print('PRQ1=', PRQ1[0])
# print('PRQ2=', PRQ2[0])
# print('PRQ4=', PRQ4[0])


label_font = 20
tick_font = 20
legend_font = 12

f_SIM_OctW = 0.96758 * 0.95  # Xmon match data
SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
f_DAC_Oct = 4.604  # 4.604
DAC_Oct_Offset = 0.5  # 0.9


Q1_base = np.mean(Q1_PSD_OctWeak[:, 1][:8])
Q2_base = np.mean(Q2_PSD_OctWeak[:, 1][:8])
Q4_base = np.mean(Q4_PSD_OctWeak[:, 1][:8])

fig, axs = plt.subplots(2, figsize=(8, 8))

f_l = 50
f_r = 535

ld1=3
ld2=5

axs_02 = axs[0].twinx()
axs_02.plot(f_SFQ, PhotonFlux, color="green", linestyle='--',
            linewidth=ld2, label='Photon Flux')
axs_02.set_ylabel("Photon flux $(\\rm m^{-2}\\rm s^{-1})$", color="black", fontsize=label_font, fontweight='bold')
axs_02.set_xlim([f_l, f_r])
axs_02.set_ylim([2e12, 6e14])
axs_02.set_yscale('log')
axs_02.legend(loc=1)
axs_02.tick_params(labelsize=tick_font)

axs[0].plot(f_SFQ, eQ1*AreaQ1, linewidth=ld2, color="red", label='$Q^{L}$')
axs[0].plot(f_SFQ, eQ4*AreaQ4, linewidth=ld2, color="black", label='$Q^{M}$')
axs[0].plot(f_SFQ, eQ2*AreaQ2, linewidth=ld2, color="blue", label='$Q^{S}$')
axs[0].set_ylabel('$e_{c}^{r}\\times A_{eff}$  $(\\rm m ^2)$', color="black", fontsize=label_font, fontweight='bold')
axs[0].set_xlim([f_l, f_r])
axs[0].set_ylim([2e-11, 6e-9])
axs[0].set_yscale('log')
axs[0].legend(loc=2)
axs[0].tick_params(labelsize=tick_font)

axs[1].axhline(y=Q1_base, color='r', linestyle='--', linewidth=ld1)
axs[1].axhline(y=Q2_base, color='b', linestyle='--', linewidth=ld1)
axs[1].axhline(y=Q4_base, color='k', linestyle='--', linewidth=ld1)

axs[1].plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', linewidth=ld2, label='$Q^{L}$')
axs[1].plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1], 'k-', linewidth=ld2, label='$Q^{M}$')
axs[1].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1], 'b-', linewidth=ld2, label='$Q^{S}$')
axs[1].tick_params(labelsize=tick_font)

axs[1].set_xlim([f_l, f_r])
axs[1].set_ylim([1e1, 1e4])
axs[1].set_yscale('log')
axs[1].legend(loc=4)
axs[1].set_xlabel("Radiator Frequency (GHz)", color="black", fontsize=label_font, fontweight='bold')
axs[1].set_ylabel("$\Gamma_{P}$ ($s^{-1}$)", color="black", fontsize=label_font, fontweight='bold')

path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
plt.tight_layout()
plt.savefig(path + '\Circmon.pdf', bbox_inches='tight', format='pdf', dpi=1200)
plt.show()