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
    X_QP = J1.Radiator["X_QP"]

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

    fig, axs = plt.subplots(2, sharex='col', figsize=(6, 4),
                            gridspec_kw={'hspace': 0.1})


    freq_l = 175
    freq_r = 310

    axs[0].plot(f, X_QP, 'k-', linewidth=4,
                label='$I_{c}$')
    axs[0].set_xlim([100, 700])
    axs[0].set_ylim([-0.02, 0.25])
    # axs[0].set_ylabel('$x_{\mathrm{QP}}$', fontsize=label_font)
    axs[0].tick_params(labelsize=tick_font)

    axs[1].plot(f, Ic_f*1e9, 'k-', linewidth=4,
                label='$I_{c}$')
    axs[1].set_xlim([150, 650])
    axs[1].set_ylim([6.5, 9.5])
    # axs[1].set_ylabel('$I_{c}$ (nA)', fontsize=label_font)
    axs[1].tick_params(labelsize=tick_font)
    axs[1].set_xlabel('Radiator frequency (GHz)', fontsize=label_font)
    axs[1].set_xticks([200, 300, 400, 500, 600])

    fig.align_ylabels(axs[:])

    path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
    plt.savefig(path+'\XmonRadiatorIcXqp.pdf', format='pdf', bbox_inches='tight', dpi=1200)
    plt.show()


