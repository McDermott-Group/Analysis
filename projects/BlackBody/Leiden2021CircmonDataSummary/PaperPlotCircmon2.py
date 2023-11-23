import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, getPhotonRate
import numpy as np
import matplotlib
#matplotlib.use("QtAgg")

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
    fileCircTest = "testpad_circmon1THz.txt"    # with wirebonds and RC and short in other ports
    fileCircTest1 = "testpad_circmon_lump-ele_wirebond.txt"    # with wirebonds and all RC in other ports
    # fileCircTest2 = "testpad_circmon_lump-ele_nowirebond.txt"    # without wirebonds, with RC in other ports
    fileQ1 = "Q1_full-chip.txt"
    fileQ2 = "Q2_full-chip.txt"
    fileQ3 = "Q3_full-chip.txt"
    fileQ4 = "Q4_full-chip.txt"

    Circ_test = AntennaCoupling()
    Circ_test.import_data(fileCircTest, JCirc, C_eff=C_eff_circ)
    f_Circ_test = Circ_test.Antenna["f"]
    ecCirc_test = Circ_test.Antenna["e_c_dB"]
    eCirc_test = Circ_test.Antenna["e_c"]
    PhotonFlux_CircTest = Circ_test.Al_Wall["PhotonFlux"]

    Circ_test1 = AntennaCoupling()
    Circ_test1.import_data(fileCircTest1, JCirc, C_eff=C_eff_circ)
    f_Circ_test1 = Circ_test1.Antenna["f"]
    ecCirc_test1 = Circ_test1.Antenna["e_c_dB"]
    eCirc_test1 = Circ_test1.Antenna["e_c"]

    # Circ_test2 = AntennaCoupling()    # with wire bonds and without, no difference
    # Circ_test2.import_data(fileCircTest2, JCirc, C_eff=C_eff_circ)
    # f_Circ_test2 = Circ_test2.Antenna["f"]
    # ecCirc_test2 = Circ_test2.Antenna["e_c_dB"]
    # eCirc_test2 = Circ_test2.Antenna["e_c"]

    Circ = AntennaCoupling()
    Circ.import_data(fileCirc, JCirc, C_eff=C_eff_circ)
    f_Circ = Circ.Antenna["f"]
    ecCirc = Circ.Antenna["e_c_dB"]
    eCirc = Circ.Antenna["e_c"]

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
    # f_SFQ = f_SFQ * f_scale / 1e9
    f_Q2 = f_Q2 * f_scale / 1e9



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


if 1:  # parity data cross talk for Aug
    f_DAC_Aug = 4.604*0.95  # 4.604
    DAC_Aug_Offset = -1
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

    fig, axs = plt.subplots(2, figsize=(8, 8), sharex=True)

    f_l = 50
    f_r = 550



    # axs_12 = axs[0].twinx()
    # axs_12.plot(f_Circ_test, PhotonFlux_CircTest, color="black", linestyle='--',
    #             linewidth=4, label='Photon Flux')
    # axs_12.set_ylabel("Photon flux $=S/hf$ $(m^{-2} sec^{-1})$", color="black", fontsize=10)
    # axs_12.set_xlim([f_l, f_r])
    # axs_12.set_ylim([5e12, 5e15])
    # axs_12.set_yscale('log')
    # axs_12.legend(loc=4)

    ld2 = 3

    # axs[0].plot(f_Circ, eCirc, color="green", linestyle='--', linewidth=ld2)
    axs[0].plot(f_Circ_test, eCirc_test, color="g", linestyle='--', linewidth=ld2)
    # axs[0].plot(f_Circ_test, eCirc_test1, color="g", marker='.', linewidth=ld2)
    axs[0].plot(f_Circ_test, eQ1, color="red", linewidth=ld2)
    axs[0].plot(f_Circ_test, eQ2, color="blue", linewidth=ld2)
    axs[0].plot(f_Circ_test, eQ4, color="black", linewidth=ld2)
    axs[0].set_xlim([f_l, f_r])
    axs[0].set_ylim([5e-5, 8e-1])
    axs[0].set_yscale('log')

    axs[1].plot((Q1_PSD_Aug[:, 0]+DAC_Aug_Offset)*f_DAC_Aug, Q1_PSD_Aug[:, 1], 'r-', linewidth=ld2)
    axs[1].plot((Q2_PSD_Aug[:, 0]+DAC_Aug_Offset)*f_DAC_Aug, Q2_PSD_Aug[:, 1], 'b-', linewidth=ld2)
    axs[1].plot((Q4_PSD_Aug[:, 0]+DAC_Aug_Offset)*f_DAC_Aug, Q4_PSD_Aug[:, 1], 'k-', linewidth=ld2)

    axs[1].set_xlim([f_l, f_r])
    axs[1].set_yscale('log')
    axs[1].set_ylim([3e1, 3e4])
    plt.show()