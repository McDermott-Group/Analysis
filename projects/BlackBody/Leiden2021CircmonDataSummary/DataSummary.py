import matplotlib.pyplot as plt
import numpy as np
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

plt.yscale('log')
plt.show()