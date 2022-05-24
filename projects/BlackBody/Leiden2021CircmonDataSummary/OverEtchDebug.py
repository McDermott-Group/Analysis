import matplotlib.pyplot as plt
from antennalib import AntennaCoupling
import numpy as np

if 1:
    """
    Import CST files and Junction parameters
    """
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff_circ = 0.5*75 * 1e-21  # Commonly used (50-100)
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
    # fileSFQ = "testpad_1.5THz.txt"
    fileSFQ_4 = "SFQ_4rectmons.txt"
    fileSFQ = "SFQ_1THz.txt"
    # fileQ_l = "Q_l_with-leads_1.5THz.txt"
    fileQ1 = "Q1_full-chip.txt"
    # fileQ2 = "Q2_full-chip.txt"
    fileQ2_overetch = "Q2_overetch.txt"
    fileQ2_non_overetch = "Q2_non_overetch.txt"
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

    # Q2 = AntennaCoupling()
    # Q2.import_data(file, JQ2, C_eff=C_eff_circ)
    # f_Q2 = Q2.Antenna["f"]
    # eQ2 = Q2.Antenna["e_c"]
    # ecQ2 = Q2.Antenna["e_c_dB"]
    # AreaQ2 = Q2.Receiver["Area"]

    Q2_oe = AntennaCoupling()
    Q2_oe.import_data(fileQ2_overetch, JQ2, C_eff=C_eff_circ)
    f_Q2_oe = Q2_oe.Antenna["f"]
    eQ2_oe = Q2_oe.Antenna["e_c"]
    ecQ2_oe = Q2_oe.Antenna["e_c_dB"]
    AreaQ2_oe = Q2_oe.Receiver["Area"]

    Q2_noe = AntennaCoupling()
    Q2_noe.import_data(fileQ2_non_overetch, JQ2, C_eff=C_eff_circ)
    f_Q2_noe = Q2_noe.Antenna["f"]
    eQ2_noe = Q2_noe.Antenna["e_c"]
    ecQ2_noe = Q2_noe.Antenna["e_c_dB"]
    AreaQ2_noe = Q2_noe.Receiver["Area"]

    Q4 = AntennaCoupling()
    Q4.import_data(fileQ4, JQ4, C_eff=C_eff_circ)
    eQ4 = Q4.Antenna["e_c"]
    ecQ4 = Q4.Antenna["e_c_dB"]
    AreaQ4 = Q4.Receiver["Area"]

    f_scale = np.sqrt(e_eff / 6.0)
    f_Circ = f_Circ * f_scale / 1e9
    f_SFQ = f_SFQ * f_scale / 1e9
    f_Q2_oe = f_Q2_oe * f_scale / 1e9
    f_Q2_noe = f_Q2_noe * f_scale / 1e9

    if 1: # receiver
        # plt.plot(f_SFQ, eQ1, c='r', label='Q1')
        plt.plot(f_Q2_oe, eQ2_oe, c='b', label='Q2 overetech')
        plt.plot(f_Q2_noe, eQ2_noe, c='k', label='Q2 non overetch')
        # plt.plot(f_SFQ, eQ4, c='k', label='Q4')
        # plt.ylim([5e-3, 8e-2])

    plt.yscale('log')
    plt.xlim([50, 1000])
    plt.ylim([1e-3, 1e-1])
    plt.legend()
    # plt.grid(which='both')
    plt.show()
