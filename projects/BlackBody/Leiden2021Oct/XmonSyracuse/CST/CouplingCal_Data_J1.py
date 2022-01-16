from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

if 1: # import CST data
    ### parameters to be tuned
    e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
    C_eff = 90 * 1e-21  # Commonly used (50-100)
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

    fileJ1 = "xmon_full-chip_JJ1.txt"
    # fileJ1 = "Xmon_testpad_with_wirebonds.txt"
    fileQ1 = "xmon_full-chip_Q1.txt"
    fileQ2 = "xmon_full-chip_Q2.txt"

    J1 = AntennaCoupling()
    J1.import_data(fileJ1, JJ1, C_eff=C_eff)
    f = J1.Antenna["f"]
    ecJ1 = J1.Antenna["e_c_dB"]
    eJ1 = J1.Antenna["e_c"]
    pgJ1 = J1.Radiator["Gamma_rad"]
    x_qpJ1 = J1.Radiator["X_QP"]
    refJ1 = J1.Al_Wall["Ref"]
    PhotonFlux = J1.Al_Wall["PhotonFlux"]

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

    f_scale = np.sqrt(e_eff / 6.0)
    f = f / 1e9
    f = f * f_scale
    f_Q1 = f_Q1 / 1e9
    f_Q1 = f_Q1 * f_scale
    f_Q2 = f_Q2 / 1e9
    f_Q2 = f_Q2 * f_scale

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
    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)
    pgJ1Q2 = []
    x_qp2photon = []
    Gamma_re = []   # photon received
    Gamma_withBase = []  # total parity rate with baseline
    l_i = 50
    for i in range(len(Area)):
        Gamma_re.append(0.5*PhotonFlux[i]*Area[i]*eQ2[i])

    ratio = 1.0/6   # for 190 GHz peak, 3.6, for 270 GHz, 7.0
    base = 110
    for i in range(len(pgJ1)):
        Gamma_withBase.append(ratio*Gamma_re[i] + base)
        # x_qp2photon.append(x_qpJ1[i]**2*5000)

    fig, axs = plt.subplots(2, 2, sharex='col', figsize=(12, 8),
                            gridspec_kw={'width_ratios': [3, 2], 'height_ratios': [2, 3],
                                         'hspace': 0.1, 'wspace': 0.05})

    axs[0, 0].plot(f, eJ1, color="red", marker="o", markersize=4,
                label='Radiator')
    axs[0, 0].plot(f_Q2, eQ2, color="blue", marker="o", markersize=4,
                label='Receiver')
    axs[0, 0].plot(f_Q2, eJ1*eQ2, color="purple", marker="o", markersize=4,
                label='Total')
    axs[0, 0].set_ylabel('Coupling Efficiency', color="black", fontsize=10)
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_xlim([l_i, 600])
    axs[0, 0].set_ylim([1e-6, 2e-1])
    axs[0, 0].legend(loc=4)
    # axs[0, 0].set_xlabel("Antenna Frequency (GHz)", color="black", fontsize=10)

    axs[1, 0].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, 'k-', label='$\Gamma_{Measured}$')
    axs[1, 0].axhline(y=base, color="blue", linestyle='--', label='$\Gamma_{0}$')
    axs[1, 0].plot(f_Q2[l_i:], Gamma_withBase[l_i:], color="grey", linestyle='-',
                label='$\Gamma_{0}+\Gamma_{PAT}$')
    # axs[1, 0].plot(f_Q2[l_i:], pgJ1Q2_scaled[l_i:],
    #             label='$photon generation rate$')
    # axs[1].plot(f_Q2[l_i:], x_qp2photon[l_i:], color="grey", linestyle='-',
    #             label='x_qp')

    axs[1, 0].set_xlabel("Radiator Josephson Frequency (GHz)", color="black",
                      fontsize=10)
    axs[1, 0].set_ylabel("$\Gamma_{p}$ ($s^{-1}$)", color="black", fontsize=10)
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlim([50, 700])
    axs[1, 0].set_ylim([100, 1100])
    axs[1, 0].legend(loc=1)
    # axs[1, 0].share

    freq_l = 175
    freq_r = 310


    axs[0, 1].plot(f, eJ1, color="red", marker="o", markersize=4,
                label='Radiator Efficiency')
    axs[0, 1].plot(f_Q2, eQ2, color="blue", marker="o", markersize=4,
                label='Receiver Efficiency')
    axs[0, 1].plot(f_Q2, eJ1 * eQ2, color="purple", marker="o", markersize=4,
                label='Total Efficiency')
    axs[0, 1].set_yscale('log')
    axs[0, 1].set_xlim([freq_l, freq_r])
    axs[0, 1].set_ylim([1e-6, 2e-1])
    axs[0, 1].yaxis.tick_right()
    # axs[0, 1].set_xlabel("Antenna Frequency (GHz)", color="black", fontsize=10)

    axs[1, 1].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, 'k-', label='$\Gamma_{Measured}$')
    # axs[1, 1].plot([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, color="black", marker="o", label='$\Gamma_{Measured}$')
    axs[1, 1].plot(f_Q2, Gamma_withBase, color="grey", linestyle='-',
                label='$\Gamma_{0}+\Gamma_{PAT}$')
    # axs[1].errorbar([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate,
    #                 yerr=Q2_J1ParityUncertainty, label='Q2', ecolor='k',
    #                 capthick=4, color='g')
    axs[1, 1].set_xlabel("Radiator Josephson Frequency (GHz)", color="black",
                      fontsize=10)
    axs[1, 1].set_yscale('log')
    axs[1, 1].set_xlim([freq_l, freq_r])
    axs[1, 1].set_ylim([100, 1100])
    axs[1, 1].yaxis.tick_right()
    # axs[1, 1].legend(loc=4)

"""
P1 Q2 and All PSD
"""
if 0:
    # plt.figure(0)
    # plt.rcParams["figure.figsize"] = (6, 20)

    fig, axs = plt.subplots(2, figsize=(7, 8))

    axs[0].errorbar(Q2_P1[:, 0] * f_AWG, Q2_P1[:, 1] * 100,
                    yerr=Q2_P1[:, 2] * 100/ np.sqrt(100), label='Q2', ecolor='g',
                    capthick=4, color='k')
    axs[0].set_xlabel("Radiator Josephson Frequency (GHz)", color="black",
                      fontsize=10)
    axs[0].set_ylabel("P1 (%)", color="black", fontsize=10)
    # axs[2].set_yscale('log')
    axs[0].set_xlim([50, 700])
    axs[0].set_ylim([1.5, 5])
    axs[0].grid(True, which="both")
    axs[0].legend(loc=4)

    axs[1].errorbar([J * f_SIM for J in Q1_J1Bias], Q1_J1ParityRate, yerr=Q1_J1ParityUncertainty, ecolor='r', capthick=4, color='g', label='$Q1$')
    axs[1].errorbar([J * f_SIM for J in Q2_J1Bias], Q2_J1ParityRate, yerr=Q2_J1ParityUncertainty, ecolor='k', capthick=4, color='k', label='$Q2$')
    axs[1].errorbar([J * f_SIM for J in Q3_J1Bias], Q3_J1ParityRate, yerr=Q3_J1ParityUncertainty, ecolor='b', capthick=4, color='b', label='$Q3$')

    axs[1].set_xlabel("Radiator Josephson Frequency (GHz)", color="black",
                      fontsize=10)
    axs[1].set_ylabel("$\Gamma_{p}$ ($s^{-1}$)", color="black", fontsize=10)
    axs[1].set_yscale('log')
    axs[1].set_xlim([50, 700])
    axs[1].set_ylim([80, 800])
    # axs[1].grid(True, which="both")
    axs[1].grid(True, which="both")
    axs[1].legend(loc=4)

plt.show()

