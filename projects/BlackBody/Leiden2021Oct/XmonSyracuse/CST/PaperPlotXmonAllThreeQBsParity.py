from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager
import os
import matplotlib
#matplotlib.use("QtAgg")

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

"""
P1 Q2 and All PSD
"""
if 1:   # one plot
    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)

    label_font = 16
    tick_font = 16
    # legend_font = 12

    # fig, axs = plt.subplots(2, sharex='col', figsize=(6, 6.5))
    plt.figure(figsize=(6, 4))

    # plt.plot(f, 100+2000 * np.multiply(x_qpJ1, x_qpJ1), 'k--', linewidth=4,
    #             label='Base+'+'$a*x_{\mathrm{qp}}^2$')

    # plt.plot(f, 100+2000 * x_qpJ1, 'k--', linewidth=4,
    #             label='Base+'+'$a*x_{\mathrm{qp}}^2$')

    plt.plot([J * f_SIM for J in Q1_J1Bias], Q1_J1ParityRate, color='r',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_1}$')
    plt.plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color='k',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_2}$')
    plt.plot([J * f_SIM for J in Q3_J1Bias], Q3_J1ParityRate, color='b',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_3}$')
    plt.yscale('log')
    plt.xlim([50, 620])
    plt.ylim([80, 1000])

    # plt.xlim([50, 1200])
    # plt.ylim([80, 3000])
    # axs[1].grid(True, which="both")
    plt.legend(loc=2, prop={'size': 15}, frameon=False)
    plt.tick_params(labelsize=tick_font)
    plt.tick_params(axis="x", direction="in", which='both')
    plt.tick_params(axis="y", direction="in", which='both')
    plt.tick_params(axis="x", width=1, length=6, which='both')
    plt.tick_params(axis="y", width=1, length=3, which='minor')
    plt.tick_params(axis="y", width=1, length=6, which='major')

    plt.xlabel("Transmitter Josephson frequency (GHz)", color="black",
                      fontsize=label_font)

    axs = plt.gca()
    axs.axvspan(368,620, facecolor='0.2', alpha=0.1)


    def fj_to_delta(f):
        return f / 92

    def delta_to_fj(delta):
        return 92 * delta

    secax = axs.secondary_xaxis('top', functions=(fj_to_delta, delta_to_fj))
    secax.tick_params(labelsize=tick_font)
    secax.set_xticks(np.arange(0, 10, 1))
    secax.set_xlabel("Transmitter Voltage Bias ($\Delta/e$)", fontsize=label_font, labelpad=10)
    secax.tick_params(axis="x", direction="in", which='both')
    secax.tick_params(axis="y", direction="in", which='both')

    secax.tick_params(axis="x", width=1, length=6, which='both')
    secax.tick_params(axis="y", width=1, length=3, which='minor')
    secax.tick_params(axis="y", width=1, length=6, which='major')

    # fig.align_ylabels(axs)
    plt.tight_layout()
    path = '..'# 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    #plt.savefig(path+'\XmonQ123Parity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.savefig('XmonQ123Parity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

if 0:   # one plot+zoom in
    label_font = 16
    tick_font = 16

    fig, axs = plt.subplots(2, figsize=(6.5, 7.5))

    axs[0].plot([J * f_SIM for J in Q1_J1Bias], Q1_J1ParityRate, color='r',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_1}$')
    axs[0].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color='k',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_2}$')
    axs[0].plot([J * f_SIM for J in Q3_J1Bias], Q3_J1ParityRate, color='b',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_3}$')
    axs[0].set_yscale('log')
    # axs[0].set_xlim([50, 620])
    # axs[0].set_ylim([80, 1000])
    axs[0].set_xlim([150, 620])
    axs[0].set_ylim([80, 1000])
    axs[0].legend(loc=2, prop={'size': 15}, frameon=False)
    axs[0].tick_params(labelsize=tick_font)
    axs[0].tick_params(axis="x", direction="in", which='both')
    axs[0].tick_params(axis="y", direction="in", which='both')
    axs[0].tick_params(axis="x", width=1, length=6, which='both')
    axs[0].tick_params(axis="y", width=1, length=3, which='minor')
    axs[0].tick_params(axis="y", width=1, length=6, which='major')

    axs[0].set_xlabel("Transmitter Josephson frequency (GHz)", color="black",
                      fontsize=label_font)

    # axs[1].plot([J * f_SIM for J in Q1_J1Bias], Q1_J1ParityRate, color='r',
    #             marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_1}$')
    # axs[1].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color='k',
    #             marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_2}$')
    axs[1].plot([J * f_SIM for J in Q3_J1Bias], Q3_J1ParityRate, color='b',
                marker="o", markersize=4, linestyle='None', label='$\mathrm{Q_3}$')
    axs[1].set_yscale('log')
    axs[1].set_xlim([245, 311])
    axs[1].set_ylim([140, 705])
    axs[1].set_xticks([250, 270, 290, 310])
    # axs[1].set_yticks([200, 300, 400, 500, 600, 700])
    # axs[1].legend(loc=2, prop={'size': 15}, frameon=False)
    axs[1].tick_params(axis="x", labelsize=tick_font)
    axs[1].tick_params(axis="y", labelsize=tick_font, which='both')
    axs[1].tick_params(axis="x", direction="in", which='both')
    axs[1].tick_params(axis="y", direction="in", which='both')
    axs[1].tick_params(axis="x", width=1, length=6, which='both')
    axs[1].tick_params(axis="y", width=1, length=3, which='minor')
    axs[1].tick_params(axis="y", width=1, length=6, which='major')

    axs[1].set_xlabel("Transmitter Josephson frequency (GHz)", color="black",
                      fontsize=label_font)




    # axs[1].set_ylabel("$\Gamma_{\mathrm{P}}$ (s$^{-1}$)", color="black", fontsize=label_font)


    # fig.align_ylabels(axs)
    plt.tight_layout()
    # path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    # plt.savefig(path+'\XmonQ123Parity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.savefig('XmonQ123Parity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()

