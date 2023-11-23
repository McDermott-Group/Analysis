from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
import os


if 1: # import CST data
    ### parameters to be tuned
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff = 75 * 1e-21  # Commonly used (50-100) 75, 88
    Jbias_offset = 1  # mDAC should be +-1 mDAC basically +-5 GHz
    k = 1  # Coupling between radiator and receiver, this could be larger than one due to the
    # fact we can generate QPs locally at the recevier's test pad

    k1 = 0.0  # coupling between on chip phonon mitigation
    f_SIM = 0.96758 * 1
    f_AWG = 4.604 * 1.05  # P1 AWG Bias freq conversion
    ### parameteres tuned done

    # JJ7 = [4.2*1e3, None, 0, 1000*150]   #[R, L, C, A] strong radiator
    JJ1 = [33 * 1e3, None, 0, 180 * 150, "Radiator"]  # [R, L, C, A] weak radiator
    JQ1 = [17.1 * 1e3, None, 0, 360 * 150, "Receiver"]
    JQ2 = [16.6 * 1e3, None, 0, 360 * 150, "Receiver"]  #
    JQ3 = [16.1 * 1e3, None, 0, 360 * 150, "Receiver"]  #

    fileJ1 = "xmon_full-chip_JJ1.txt"
    # fileJ1 = "Xmon_testpad_with_wirebonds.txt"
    fileQ1 = "xmon_full-chip_Q1.txt"
    fileQ2 = "xmon_full-chip_Q2.txt"
    fileQ3 = "xmon_full-chip_Q3.txt"

    J1 = AntennaCoupling()
    J1.import_data(fileJ1, JJ1, C_eff=C_eff)
    f = J1.Antenna["f"]
    ecJ1 = J1.Antenna["e_c_dB"]
    eJ1 = J1.Antenna["e_c"]
    pgJ1 = J1.Radiator["Gamma_rad"]
    x_qpJ1 = J1.Radiator["X_QP"]
    refJ1 = J1.Al_Wall["Ref"]
    PhotonFlux = J1.Al_Wall["PhotonFlux"]
    Ic_f = J1.Radiator["Ic_f"]

    Q1 = AntennaCoupling()
    Q1.import_data(fileQ1, JQ1, C_eff=C_eff)
    f_Q1 = Q1.Antenna["f"]
    eQ1 = Q1.Antenna["e_c"]

    Q2 = AntennaCoupling()
    Q2.import_data(fileQ2, JQ2, C_eff=C_eff)
    f_Q2 = Q2.Antenna["f"]
    ecQ2 = Q2.Antenna["e_c_dB"]
    eQ2 = Q2.Antenna["e_c"]
    Area = Q2.Receiver["Area"]

    Q3 = AntennaCoupling()
    Q3.import_data(fileQ3, JQ3, C_eff=C_eff)
    f_Q3 = Q3.Antenna["f"]
    eQ3 = Q3.Antenna["e_c"]

    f_scale = np.sqrt(e_eff / 6.0)
    f = f / 1e9
    f = f * f_scale
    f_Q1 = f_Q1 / 1e9
    f_Q1 = f_Q1 * f_scale
    f_Q2 = f_Q2 / 1e9
    f_Q2 = f_Q2 * f_scale
    f_Q3 = f_Q3 / 1e9
    f_Q3 = f_Q3 * f_scale

    # plt.plot(f, eQ1, color="red", marker="o", markersize=4, label='Q1')
    # plt.plot(f, eQ2, color="blue", marker="o", markersize=4, label='Q2')
    # plt.plot(f, eQ3, color="black", marker="o", markersize=4, label='Q3')
    #
    # plt.legend()
    # plt.show()

if 1:
    """
    Import measurement data starts
    """
    Q1_PSD_Data_file = "Q1_PSD_Data.txt"
    Q2_PSD_Data_file = "Q2_PSD_Data.txt"
    Q3_PSD_Data_file = "Q3_PSD_Data.txt"

    Q2_P1_Data_file = "Q2_P1_Data.txt"

    Q1_J1Bias = np.loadtxt(Q1_PSD_Data_file, usecols=[0], skiprows=3)
    Q1_J1ParityRate = np.loadtxt(Q1_PSD_Data_file, usecols=[1], skiprows=3)
    Q1_J1ParityUncertainty = np.loadtxt(Q1_PSD_Data_file, usecols=[2], skiprows=3)
    Q1_J1Fidelity = np.loadtxt(Q1_PSD_Data_file, usecols=[3], skiprows=3)

    Q2_J1Bias = np.loadtxt(Q2_PSD_Data_file, usecols=[0], skiprows=3)
    Q2_J1ParityRate = np.loadtxt(Q2_PSD_Data_file, usecols=[1], skiprows=3)
    Q2_J1ParityUncertainty = np.loadtxt(Q2_PSD_Data_file, usecols=[2], skiprows=3)
    Q2_J1Fidelity = np.loadtxt(Q2_PSD_Data_file, usecols=[3], skiprows=3)

    Q3_J1Bias = np.loadtxt(Q3_PSD_Data_file, usecols=[0], skiprows=3)
    Q3_J1ParityRate = np.loadtxt(Q3_PSD_Data_file, usecols=[1], skiprows=3)
    Q3_J1ParityUncertainty = np.loadtxt(Q3_PSD_Data_file, usecols=[2], skiprows=3)
    Q3_J1Fidelity = np.loadtxt(Q3_PSD_Data_file, usecols=[3], skiprows=3)

    Q2_P1 = np.loadtxt(Q2_P1_Data_file, usecols=[0, 1, 2], skiprows=3)

    """
    Import measurement data ends
    """


"""
Calculate the noise bandwidth
"""

"""
Q2 # polished
"""
if 1:

    label_font = 16
    tick_font = 13
    legend_font = 15

    legend_font = font_manager.FontProperties(
        # weight='bold',
        style='normal', size=legend_font)

    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)
    pgJ1Q2 = []
    x_qp2photon = []
    Gamma_re = []   # photon received
    Gamma_withBase = []  # total parity rate with baseline
    l_i = 50
    r_i = 700
    for i in range(len(Area)):
        Gamma_absorbed = 0.5*PhotonFlux[i]*Area[i]*eQ2[i]
        Gamma_re.append(Gamma_absorbed)

    ratio = 0.5   # for 190 GHz peak, 3.6, for 270 GHz, 7.0
    base = 110
    for i in range(len(pgJ1)):
        if f[i] >= 92:
            Gamma_withBase.append(ratio*Gamma_re[i] + base)
        else:
            Gamma_withBase.append(base)

    fig, axs = plt.subplots(2, 2, figsize=(12, 8),
                            gridspec_kw={'width_ratios': [3, 2], 'height_ratios': [2, 3],
                                         'hspace': 0.1, 'wspace': 0.1})

    axs[0, 0].plot(f, eJ1, color="red", marker="o", markersize=4,
                label='Transmitter')
    axs[0, 0].plot(f_Q2, eQ2, color="blue", marker="o", markersize=4,
                label='Receiver')
    axs[0, 0].plot(f_Q2, eJ1*eQ2, color="purple", marker="o", markersize=4,
                label='Total')
    axs[0, 0].set_ylabel('$e_{\mathrm{c}}$', color="black", fontsize=label_font)
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_xlim([l_i, r_i])
    axs[0, 0].set_ylim([1e-6, 2e-1])
    axs[0, 0].tick_params(labelsize=tick_font)
    axs[0, 0].legend(loc=4, prop=legend_font)

    axs[1, 0].plot(f_Q2, Gamma_withBase, color="grey", linestyle='-',  linewidth=2,
                label='$\Gamma_{0}+\Gamma_{\mathrm{PAT}}$')
    axs[1, 0].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color="black", linewidth=2,
                   marker="o", markersize=6, label='$\Gamma_{\mathrm{tot}}$')
    # axs[1, 0].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color="black", linestyle='d', markersize=4, label='$\Gamma_{Measured}$')
    axs[1, 0].axhline(y=base, c='g', linestyle='--', label='$\Gamma_{0}$')

    axs[1, 0].set_xlabel("Radiator Frequency (GHz)", color="black",
                      fontsize=label_font)
    axs[1, 0].set_ylabel("$\Gamma_{P}$ ($s^{-1}$)", color="black", fontsize=label_font)
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlim([l_i, r_i])
    axs[1, 0].set_ylim([100, 1100])
    axs[1, 0].legend(loc=1, prop=legend_font, frameon=False)
    axs[1, 0].tick_params(labelsize=tick_font)
    # axs[1, 0].share

    freq_l = 175
    freq_r = 310


    axs[0, 1].plot(f, Ic_f*1e9, color="black", marker="o", markersize=4,
                label='$I_{c}$')
    axs[0, 1].set_xlim([100, 700])
    axs[0, 1].set_ylim([6.5, 9.5])
    axs[0, 1].set_ylabel('$I_{c}$ (nA)', fontsize=label_font)
    axs[0, 1].yaxis.set_label_position("right")
    axs[0, 1].tick_params(labelsize=tick_font)
    axs[0, 1].yaxis.tick_right()

    axs[1, 1].plot(f_Q2, Gamma_withBase, color="grey", linestyle='-',  linewidth=4,
                label='$\Gamma_{0}+\Gamma_{PAT}$')
    axs[1, 1].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color="black", linewidth=4,
                   marker="o", markersize=6, label='$\Gamma_{Measured}$')
    axs[1, 1].set_xlabel("Radiator Frequency (GHz)", color="black",
                      fontsize=label_font)
    axs[1, 1].set_yscale('log')
    axs[1, 1].set_xlim([freq_l, freq_r])
    axs[1, 1].set_ylim([100, 1100])
    axs[1, 1].tick_params(labelsize=tick_font)
    axs[1, 1].yaxis.tick_right()
    axs[1, 1].set_ylabel("$\Gamma_{P}$ ($s^{-1}$)", color="black", fontsize=label_font)
    axs[1, 1].yaxis.set_label_position("right")

    fig.align_ylabels(axs)

    # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    # plt.savefig(path+'\XmonSpectroscopy.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

"""
P1 Q2 and All PSD
"""
if 1:
    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)

    label_font = 16
    tick_font = 13
    # legend_font = 12

    fig, axs = plt.subplots(2, sharex='col', figsize=(6, 6.5))

    axs[0].plot([J * f_SIM for J in Q1_J1Bias], Q1_J1ParityRate, color='r',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_1}$')
    axs[0].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color='k',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_2}$')
    axs[0].plot([J * f_SIM for J in Q3_J1Bias], Q3_J1ParityRate, color='b',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_3}$')
    axs[0].set_yscale('log')
    axs[0].set_xlim([50, 620])
    axs[0].set_ylim([80, 1000])
    # axs[1].grid(True, which="both")
    axs[0].legend(loc=2, prop={'size': 15}, frameon=False)
    axs[0].tick_params(labelsize=tick_font)
    axs[0].tick_params(axis="x", direction="in", which='both')
    axs[0].tick_params(axis="y", direction="in", which='both')
    axs[0].tick_params(axis="x", width=1, length=6, which='both')
    axs[0].tick_params(axis="y", width=1, length=3, which='minor')
    axs[0].tick_params(axis="y", width=1, length=6, which='major')

    axs[1].errorbar(Q2_P1[:, 0] * f_AWG, Q2_P1[:, 1] *100, yerr=Q2_P1[:, 2]*100 / np.sqrt(100), label='$Q_{2}$', ecolor='k',
                    capthick=4, color='k', linewidth=2, fmt="o")
    # axs[0].set_ylabel("$P_{1}$ (%)", color="black", fontsize=label_font)
    # axs[2].set_yscale('log')
    axs[1].set_xlim([50, 620])
    axs[1].set_ylim([1.8, 4.5])
    axs[1].tick_params(labelsize=tick_font)
    axs[1].tick_params(axis="x", direction="in", which='both')
    axs[1].tick_params(axis="y", direction="in", which='both')
    axs[1].tick_params(axis="x", width=1, length=6, which='both')
    axs[1].tick_params(axis="y", width=1, length=6, which='both')
    axs[1].set_xlabel("Transmitter frequency (GHz)", color="black",
                      fontsize=label_font)




    # axs[1].set_ylabel("$\Gamma_{\mathrm{P}}$ (s$^{-1}$)", color="black", fontsize=label_font)


    fig.align_ylabels(axs)
    plt.tight_layout()
    #path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    #plt.savefig('XmonQ123Parity_P1.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

