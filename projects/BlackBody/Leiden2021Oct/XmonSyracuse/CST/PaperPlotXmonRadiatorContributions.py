from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.font_manager as font_manager
from matplotlib.ticker import ScalarFormatter, NullFormatter, FormatStrFormatter
import os
import matplotlib
#matplotlib.use("tkAgg")

# plot for main text or for supplement?
supplement = True

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

    JJ7 = [4.2*1e3, None, 0, 1000*150, "Radiator"]   #[R, L, C, A] strong radiator
    JJ1 = [33 * 1e3, None, 0, 180 * 150, "Radiator"]  # [R, L, C, A] weak radiator
    JQ1 = [17.1 * 1e3, None, 0, 360 * 150, "Receiver"]
    # JQ2 = [16.6 * 1e3, None, 0, 360 * 150, "Receiver"]  #
    JQ2 = [16.6 * 1e3, None, 0, 390 * 156, "Receiver"]  #
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
    ShotNoiseFlux = J1.Al_Wall["PhotonFlux_ShotNoise"]
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

    label_font = 44
    tick_font = 44
    legend_font = 15

    legend_font = font_manager.FontProperties(
        # weight='bold',
        style='normal', size=legend_font)

    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)
    pgJ1Q2 = []
    x_qp2photon = []
    powers_radiated_coherent = []   # photon received
    powers_radiated_incoherent = []
    l_i = 50
    r_i = 610 if not supplement else 1000
    polarization_factor = 0.5

    for i in range(len(Area)):
        if f[i] < 92:
            powers_radiated_coherent.append(0)
            powers_radiated_incoherent.append(0)
            continue

        # coherent radiation:
        power_radiated_coherent = PhotonFlux[i]*f[i]
        # radiation from shot noise:
        power_radiated_incoherent = 0
        for j in range(len(Area)):
            power_radiated_incoherent+=ShotNoiseFlux[i][j]*f[j]

        powers_radiated_coherent.append(power_radiated_coherent*6.63e-34*1e18)
        powers_radiated_incoherent.append(power_radiated_incoherent*6.63e-34*1e18)

    # fig, axs = plt.subplots(2, figsize=(30, 22),
    #                         gridspec_kw={'height_ratios': [2, 3], 'hspace': 0.75})

    fig = plt.figure(figsize=(30, 12))#(2, figsize=(30, 12),
                            #gridspec_kw={'height_ratios': [1,0.01], 'hspace': 0.0})
    axs = fig.gca()
    #
    # axs[0].plot(f, eJ1, color="red", marker="o", markersize=4,linewidth=8,
    #             label='                ')
    # axs[0].plot(f_Q2, eQ2, color="blue", marker="o", markersize=4,linewidth=8,
    #             label='                ')
    # axs[0].plot(f_Q2, eJ1*eQ2, color="purple", marker="o", markersize=4,linewidth=8,
    #             label='                ')
    # axs[0].set_ylabel('   ', color="black", fontsize=label_font)
    # axs[0].set_yscale('log')
    # axs[0].set_xlim([l_i, r_i])
    # axs[0].set_ylim([1e-6, 2e-1])
    # axs[0].set_yticks([1e-5, 1e-3, 1e-1])
    # axs[0].tick_params(labelsize=tick_font)
    # # axs[0].legend(loc=4, prop=legend_font)
    # axs[0].tick_params(axis="x", direction="in", which='both')
    # axs[0].tick_params(axis="y", direction="in", which='both')
    # axs[0].tick_params(axis="x", width=1, length=6, which='both')
    # axs[0].tick_params(axis="y", width=1, length=6, which='major')
    # axs[0].tick_params(axis="y", width=1, length=3, which='minor')

    # axs[0].set_ylabel("$e_{\mathrm{c}}$)", color="black", fontsize=label_font)

    axs.plot(f_Q2, powers_radiated_coherent, color="red", linestyle='-',  linewidth=8)
    axs.plot(f_Q2, powers_radiated_incoherent, color="blue", linestyle='-',  linewidth=8)

    # axs[1].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, 'o', color="black")

    # axs[1].set_xlabel("Transmitter Frequency (GHz)", color="black",
    #                   fontsize=label_font)
    # axs[1].set_ylabel("$\Gamma_{\mathrm{p}}$ ($\mathrm{s^{-1}}$)", color="black", fontsize=label_font)
    axs.set_yscale('log')
    axs.set_xlim([0, 1000])
    axs.set_ylim([0.02, 200])
    # axs[1].legend(loc=2, prop=legend_font, frameon=False)
    axs.set_ylabel("Radiated Power (aW)", color="black", fontsize=label_font)
    axs.set_xlabel("Transmitter Josephson Frequency (GHz)", color="black", fontsize=label_font)
    axs.tick_params(labelsize=tick_font)

    def fj_to_delta(f):
        return f / 46 / 2

    def delta_to_fj(delta):
        return 46 * delta * 2

    axs.axvspan(368,1100, facecolor='0.2', alpha=0.1)
    secax = axs.secondary_xaxis('top', functions=(fj_to_delta, delta_to_fj))
    secax.tick_params(labelsize=tick_font)
    secax.xaxis.set_major_formatter(FormatStrFormatter('%2d'))
    secax.set_xticks(np.arange(0,27,2))
    secax.set_xlabel("Transmitter Voltage Bias ($\Delta/e$)", fontsize=label_font, labelpad=15)
    secax.tick_params(axis="x", direction="in", which='both')
    secax.tick_params(axis="y", direction="in", which='both')

    secax.tick_params(axis="x", width=1, length=6, which='both')
    secax.tick_params(axis="y", width=1, length=3, which='minor')
    secax.tick_params(axis="y", width=1, length=6, which='major')


    axs.tick_params(axis="x", direction="in", which='both')
    axs.tick_params(axis="y", direction="in", which='both')

    axs.tick_params(axis="x", width=1, length=6, which='both')
    axs.tick_params(axis="y", width=1, length=3, which='minor')
    axs.tick_params(axis="y", width=1, length=6, which='major')

    plt.savefig("XMonPowerContributions.png", format="png")
    plt.show()
