import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, getPhotonRate
import numpy as np
import matplotlib
matplotlib.use("QtAgg")

if 1:
    """
    Import CST files and Junction parameters
    """
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff_circ = 75 * 1e-21  # Commonly used (50-100)
    C_eff_SFQ = 75 * 1e-21  # Commonly used (50-100)

    JCirc = [4.8 * 1e3, None, 0, 320 * 123 * 2, "Radiator"]  # [R, L, C, A]
    JSFQ_weak = [16 * 1e3, None, 0, 100 * 200, "Radiator"]  # [R, L, C, A]
    JSFQ_weak_10Rn = [10*16 * 1e3, None, 0, 100 * 200, "Radiator"]  # [R, L, C, A]
    JSFQ_weak_100Rn = [100*16 * 1e3, None, 0, 100 * 200, "Radiator"]  # [R, L, C, A]
    JSFQ_weak_qp = [16 * 1e3 * 100, None, 0, 100 * 200, "Radiator"]  # [R, L, C, A]
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

    #Radiator:
    SFQ_4 = AntennaCoupling()
    SFQ_4.import_data(fileSFQ_4, JSFQ_weak, C_eff=C_eff_SFQ)
    ecSFQ_4 = SFQ_4.Antenna["e_c_dB"]
    eSFQ_4 = SFQ_4.Antenna["e_c"]
    PhotonFlux = SFQ_4.Al_Wall["PhotonFlux"]

    SFQ_4_10Rn = AntennaCoupling()
    SFQ_4.import_data(fileSFQ_4, JSFQ_weak_10Rn, C_eff=C_eff_SFQ)
    ecSFQ_4_10Rn = SFQ_4.Antenna["e_c_dB"]
    eSFQ_4_10Rn = SFQ_4.Antenna["e_c"]
    PhotonFlux = SFQ_4.Al_Wall["PhotonFlux"]

    SFQ_4_100Rn = AntennaCoupling()
    SFQ_4.import_data(fileSFQ_4, JSFQ_weak_100Rn, C_eff=C_eff_SFQ)
    ecSFQ_4_100Rn = SFQ_4.Antenna["e_c_dB"]
    eSFQ_4_100Rn = SFQ_4.Antenna["e_c"]
    PhotonFlux = SFQ_4.Al_Wall["PhotonFlux"]

    SFQ_4_qp = AntennaCoupling()
    SFQ_4_qp.import_data(fileSFQ_4, JSFQ_weak_qp, C_eff=C_eff_SFQ)
    ecSFQ_4_qp= SFQ_4_qp.Antenna["e_c_dB"]
    eSFQ_4_qp = SFQ_4_qp.Antenna["e_c"]
    PhotonFlux_qp = SFQ_4_qp.Al_Wall["PhotonFlux"]

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


label_font = 14
tick_font = 24
legend_font = 12

f_SIM_OctW = 0.96758 * 0.95  # Xmon match data
SIM_OctW_Offset = 2  # Xmon match data   2mV is OK or max
f_DAC_Oct = 4.604  # 4.604
DAC_Oct_Offset = 0.5  # 0.9


Q1_base = np.mean(Q1_PSD_OctWeak[:, 1][:8])
Q2_base = np.mean(Q2_PSD_OctWeak[:, 1][:8])
Q4_base = np.mean(Q4_PSD_OctWeak[:, 1][:8])

fig, axs = plt.subplots(2, figsize=(12, 14),
                        gridspec_kw={'height_ratios': [2.5, 3],
                                     'hspace': 0.15})

f_l = 50
f_r = 535

ld1=6
ld2=6


# axs[0].plot(f_SFQ, eSFQ_4, color="green", linestyle='--', linewidth=ld2, label='Transmitter')
# axs[0].plot(f_SFQ, eSFQ_4, color=[0, 1, 0], linestyle='--', linewidth=ld2, label='Transmitter')
axs[0].plot(f_SFQ, eSFQ_4, color="darkgreen", linestyle='--', linewidth=ld2, label='            ')
axs[0].plot(f_SFQ[0:185], list(eSFQ_4_10Rn[0:184])+list([eSFQ_4[185]]), color="darkgreen", linestyle='-.', linewidth=ld2, label='            ')
axs[0].plot(f_SFQ[0:185], list(eSFQ_4_100Rn[0:184])+list([eSFQ_4_10Rn[185]]), color="darkgreen", linestyle='dotted', linewidth=ld2, label='            ')

# axs[0].plot(f_SFQ, eSFQ_4_qp, color="green", linestyle='--', linewidth=ld2, label='            ')

axs[0].plot(f_SFQ, eQ1, linewidth=ld2, color="red", label='            ')
axs[0].plot(f_SFQ, eQ4, linewidth=ld2, color="black", label='   ')
axs[0].plot(f_SFQ, eQ2, linewidth=ld2, color="blue", label='   ')
axs[0].set_xlim([f_l, f_r])
axs[0].set_ylim([4e-5, 1e-1])
axs[0].set_yscale('log')
# axs[0].legend(loc=4, ncol=2, frameon=False, prop={'size': 13})
axs[0].tick_params(labelsize=tick_font)

axs[0].tick_params(axis="x", direction="in", which='both')
axs[0].tick_params(axis="y", direction="in", which='both')
axs[0].tick_params(axis="x", width=1, length=4, which='both')
axs[0].tick_params(axis="y", width=1, length=3, which='minor')
axs[0].tick_params(axis="y", width=1, length=6, which='major')

axs[1].axhline(y=Q1_base, color='r', linestyle='--', linewidth=ld1)
axs[1].axhline(y=Q2_base, color='b', linestyle='--', linewidth=ld1)
axs[1].axhline(y=Q4_base, color='k', linestyle='--', linewidth=ld1)

axs[1].plot((Q1_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q1_PSD_OctWeak[:, 1], 'r-', linewidth=ld2, label='$\mathrm{Q_{1}}$')
axs[1].plot((Q4_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q4_PSD_OctWeak[:, 1], 'k-', linewidth=ld2, label='$\mathrm{Q_{2}}$')
axs[1].plot((Q2_PSD_OctWeak[:, 0]+SIM_OctW_Offset)*f_SIM_OctW, Q2_PSD_OctWeak[:, 1], 'b-', linewidth=ld2, label='$\mathrm{Q_{3}}$')
axs[1].tick_params(labelsize=tick_font)

axs[1].set_xlim([f_l, f_r])
axs[1].set_ylim([9e0, 1e4])
axs[1].set_yscale('log')
# axs[1].legend(loc=4, ncol=3, frameon=False)
# axs[1].set_xlabel("Radiator Frequency (GHz)", color="black", fontsize=label_font)

axs[1].tick_params(axis="x", direction="in", which='both')
axs[1].tick_params(axis="y", direction="in", which='both')
axs[1].tick_params(axis="x", width=1, length=4, which='both')
axs[1].tick_params(axis="y", width=1, length=3, which='minor')
axs[1].tick_params(axis="y", width=1, length=6, which='major')

plt.tight_layout()
path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
plt.savefig(path + '\CircmonResistance.pdf', bbox_inches='tight', format='pdf', dpi=1200, transparent=True)
plt.show()