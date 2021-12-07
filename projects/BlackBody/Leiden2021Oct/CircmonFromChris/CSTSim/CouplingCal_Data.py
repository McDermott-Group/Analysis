from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np

### parameters to be tuned
e_eff = 6 # limit (1, 6.5), the voltage can also be built in to have a larger range
C_eff = 180*1e-21   # Commonly used (50-100)
Jbias_offset = 0   # mDAC should be +-1 mDAC basically +-5 GHz
k = 1   # Coupling between radiator and receiver, this could be larger than one due to the
        # fact we can generate QPs locally at the recevier's test pad

k1 = 0.0   # coupling between on chip phonon mitigation
# f_SIM = 0.97
f_SIM = 0.968
### parameteres tuned done

JSFQ = [16*1e3, None, 0, 100*200]   #[R, L, C, A]
JQ1 = [16.6*1e3, 18.3*1e-9, 0, 193.8*121.8]   #
JQ2 = [13.2*1e3, 14.6*1e-9, 0, 350*126.6]   #
JQ3 = [19*1e3, 21*1e-9, 0, 310*126]   #    some issue with Q3
JQ4 = [15*1e3, 19.9*1e-9, 0, 184.4*122.5*2]   #    some issue with Q3

fileSFQ = "testpad_1.5THz.txt"
# fileQ1 = "Q1.txt"
# fileQ1 = "Q1_with-leads_1.5THz.txt"
# fileQ2 = "Q2_1.5THz.txt"
# fileQ3 = "Q3_with-leads_1.5THz.txt"
# fileQ4 = "Q4_1.5THz.txt"

fileQ1 = "Q1_full-chip.txt"
fileQ2 = "Q2_full-chip.txt"
fileQ3 = "Q3_full-chip.txt"
fileQ4 = "Q4_full-chip.txt"

SFQ = AntennaCoupling()
SFQ.import_data(fileSFQ, JSFQ, C_eff=100*1e-21)
f_SFQ = SFQ.Antenna["f"]
ecSFQ = SFQ.e_c_dB
eSFQ = SFQ.e_c
pgSFQ = SFQ.p_g
refSFQ = SFQ.ref

# SFQ.plot()

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1, C_eff=C_eff)
ecQ1 = Q1.e_c_dB
eQ1 = Q1.e_c
Z_ReQ1 = Q1.Antenna["Z_Re"]
f_Q1 = Q1.Antenna["f"]

Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2, C_eff=C_eff)
ecQ2 = Q2.e_c_dB
eQ2 = Q2.e_c
Z_ReQ2 = Q2.Antenna["Z_Re"]

Q3 = AntennaCoupling()
Q3.import_data(fileQ3, JQ3, C_eff=C_eff)
ecQ3 = Q3.e_c_dB
eQ3 = Q3.e_c
Z_ReQ3 = Q3.Antenna["Z_Re"]

Q4 = AntennaCoupling()
Q4.import_data(fileQ4, JQ4, C_eff=C_eff)
ecQ4 = Q4.e_c_dB
eQ4 = Q4.e_c
Z_ReQ4 = Q4.Antenna["Z_Re"]


Q1_PSD = np.array([
    [0, 1009.34], [10, 958.58], [20, 1145.83], [30, 1046.75], [40, 1160.41],
    [50, 944.89], [60, 1210.19], [65, 1009.63], [70, 962.35], [75, 1140.02],
    [80, 1023.09], [85, 950.27], [90, 1115.45], [95, 1010.86], [100, 1010.15],
    [105, 941.78], [110, 1220.46], [115, 3437.57], [120, 4664.58], [125, 5172.3],
    [130, 4358.9], [135, 4078.99], [140, 3504.63], [145, 2565.52], [150, 2402.98],
    [155, 2252.7], [160, 2143.96], [165, 2867.33], [170, 2797.36], [175, 2298.0],
    [180, 2121.93], [185, 1866.51], [190, 2393.17], [195, 3000.95], [200, 3228.92],
    [205, 2731.99], [210, 4102.09], [215, 2740.09], [220, 2366.99], [225, 2074.19],
    [230, 2532.39], [235, 2125.4], [240, 1961.89], [245, 1994.41], [250, 1631.03],
    [255, 1733.2], [260, 1552.79], [265, 1505.57], [270, 1453.08], [275, 1596.24],
    [280, 1432.32], [285, 1493.91], [290,1538.86], [295, 1804.63], [300, 1614.51],
    [305, 1580.93], [310, 1528.13], [315, 1682.05], [320, 1507.64], [325, 1496.16],
    [330, 1612.51], [335, 1797.43], [340, 1642.42], [345, 1414.63], [350, 1333.24],
    [355, 1471.76], [360, 1455.64], [365, 1274.01], [370, 1213.63], [375, 1250.93],
    [380, 1262.28], [385, 1188.31], [390, 1350.99], [395, 1381.67], [400, 1434.28],
    [405, 1381.46], [410, 1514.37], [415, 1276.58], [420, 1198.76], [425, 1342.48],
    [430, 1600.78], [435, 1932.09], [440, 3354.68], [445, 2190.77], [450, 4330.97],
    [455, 4096.65], [455, 4096.65], [460, 4260.12], [465, 4579.53], [470, 5555.07],
    [475, 6303.84], [480, 5552.77], [485, 4612.55], [490, 5192.22], [495, 6452.49],
    [500, 6584.32], [505, 5684.13], [510, 5332.56], [515, 5229.72], [520, 5744.56],
    [525, 5112.43], [530, 6392.96], [535, 6664.02], [540, 5807.01], [545, 8712.34],
    [550, 5562.97], [555, 6059.17], [560, 5850.5], [565, 6664.39], [570, 7159.78],
    [575, 7547.52], [580, 7821.21]
])
Q2_PSD = np.array([
    [0, 12.86], [10, 11.36], [20, 13.26], [30, 12.21], [40, 12.7], [50, 12.83],
    [60, 12.84], [65, 11.72], [70, 12.77], [75, 13.46], [80,14.29], [85, 11.23],
    [90, 12.81], [95, 12.3], [100, 12.55], [105, 14.13], [110, 13.12], [115, 13.11],
    [120, 13.91], [125, 46.08], [130, 70.16], [135, 59.33], [140, 49.56], [145, 45.13],
    [150, 33.69], [155, 32.47], [160, 29.91], [165, 40.26], [170, 49.44], [175, 74.12],
    [180, 69.88], [185, 53.25], [190, 76.92], [195, 87.42], [200, 104.56], [205, 138.78],
    [210, 135.9], [215, 103.5], [220, 133.22], [225, 111.34], [230, 101.3], [235, 114.74],
    [240, 98.55], [245, 92.63], [250, 79.81], [255, 79.2], [260, 93.96], [265, 72.88],
    [270, 76.79], [275, 89.38], [280, 84.88], [285, 89.87], [290, 80.34], [295, 113.91],
    [300, 91.7], [305, 103.52], [310, 137.55], [315, 242.21], [320, 289.37], [325, 309.71],
    [330, 394.06], [335, 445.71], [340, 398.16], [345, 328.41], [350, 332.01], [355, 378.13],
    [360, 397.16], [365, 407.39], [370, 417.23], [375, 501.99], [380, 589.45], [385, 539.99],
    [390, 591.56], [395, 581.89], [400, 598.1], [405, 587.99], [410, 548.04], [415, 464.87],
    [420, 344.23], [425, 326.42], [430, 332.82], [435, 317.92], [440, 269.24],
    [445, 235.44], [450, 277.23], [455, 229.28], [455, 229.28], [460, 189.52],
    [465, 180.09], [470, 176.76], [475, 144.67], [480, 142.35], [485, 141.16],
    [490, 139.55], [495, 154.04], [500, 166.06], [505, 161.47], [510, 156.82],
    [515, 140.3], [520, 123.91], [525, 125.29], [530, 127.77], [535, 122.87],
    [540, 128.8], [545, 118.93], [550, 126.27], [555, 128.9], [560, 135.61],
    [565, 144.21], [570, 151.43], [575, 153.62], [580, 156.3]
])
Q4_PSD = np.array([
    [0, 201.33], [10, 191.95], [20, 185.84], [30, 184.72], [40, 190.91], [50, 192.32],
    [60, 191.89], [65, 199.26], [70, 185.41], [75, 186.42], [80, 191.55], [85, 187.11],
    [90, 196.13], [95, 192.59], [100, 183.86], [105, 195.27], [110, 173.7], [115, 197.69],
    [120, 230.62], [125, 303.48], [130, 325.11], [135, 437.51], [140, 443.43], [145, 591.88],
    [150, 542.91], [155, 349.66], [160, 333.63], [165, 492.91], [170, 650.01], [175, 874.25],
    [180, 789.24], [185, 743.13], [190, 1292.96], [195, 2104.6], [200, 2590.91],
    [205, 2378.21], [210, 3400.27], [215, 2979.89], [220, 2583.41], [225, 3639.7],
    [230, 2984.22], [235, 3062.72], [240, 2624.84], [245, 2456.55], [250, 2150.58],
    [255, 2201.17], [260, 2371.47], [265, 2309.56], [270, 2076.11], [275, 2427.54],
    [280, 2452.78], [285, 1863.37], [290, 1673.27], [295, 1555.0], [300, 1159.26],
    [305, 1230.32], [310, 1226.15], [315, 1300.76], [320, 1719.83], [325, 1759.15],
    [330, 1350.34], [335, 1180.01], [340, 1433.55], [345, 1233.66], [350, 1099.69],
    [355, 826.48], [360, 706.27], [365, 750.64], [370, 655.02], [375, 760.93],
    [380, 798.0], [385, 684.65], [390, 619.28], [395, 672.03], [400, 646.95],
    [405, 650.66], [410, 591.86], [415, 603.73], [420, 661.36], [425, 708.76],
    [430, 675.09], [435, 863.73], [440, 834.24], [445, 846.6], [450, 891.59],
    [455, 900.44], [455, 900.44], [460, 936.8], [465, 840.21], [470, 927.87],
    [475, 862.11], [480, 860.69], [485, 919.19], [490, 918.64], [495, 929.98],
    [500, 944.81], [505, 994.55], [510, 984.28], [515, 974.87], [520, 1126.25],
    [525, 1080.68], [530, 1000.69], [535, 1014.04], [540, 1088.76], [545, 1032.05],
    [550, 1028.85], [555, 1111.87], [560, 1200.51], [565, 1246.12], [570, 1264.38],
    [575, 1305.01], [580, 1357.58]
])


f_scale = np.sqrt(e_eff/6.0)
f_SFQ = f_SFQ/1e9
f_SFQ = f_SFQ*f_scale
f_Q1 = f_Q1/1e9
f_Q1 = f_Q1*f_scale

### two vertical plots
fig, axs = plt.subplots(4)
axs[0].plot(f_SFQ, ecSFQ, color="red", marker="o")
axs[0].set_ylabel("Radiator (dB)",color="black", fontsize=10)
axs[0].set_xlim([50, 600])
axs[0].set_ylim([-60, 0])
axs[0].grid()

axs[1].plot(f_Q1, ecQ1, color="blue", marker="o", label='Q1')
axs[1].plot(f_Q1, ecQ2, color="red", marker="o", label='Q2')
axs[1].plot(f_Q1, ecQ4, color="green", marker="o", label='Q4')
axs[1].set_ylabel("Receiver QB (dB)",color="black", fontsize=10)
axs[1].set_xlim([50, 600])
axs[1].grid()
axs[1].legend()
axs[1].set_ylim([-40, -10])

axs[2].plot(f_Q1, ecQ1+ecSFQ*k, color="blue", marker="o", label='Q1')
axs[2].plot(f_Q1, ecQ2+ecSFQ*k, color="red", marker="o", label='Q2')
axs[2].plot(f_Q1, ecQ4+ecSFQ*k, color="green", marker="o", label='Q4')
axs[2].set_ylabel("Receiver QB (dB)", color="black", fontsize=10)
axs[2].set_xlim([50, 600])
axs[2].grid()
axs[2].legend()
axs[2].set_ylim([-60, -20])


axs[3].plot(Q1_PSD[:, 0]*f_SIM, Q1_PSD[:, 1], color='b', label='Q1_PSD')
axs[3].plot(Q2_PSD[:, 0]*f_SIM, Q2_PSD[:, 1], color='r', label='Q1_PSD')
axs[3].plot(Q4_PSD[:, 0]*f_SIM, Q4_PSD[:, 1], color='g', label='Q1_PSD')
axs[3].set_xlabel("Freq (GHz)", color="black", fontsize=10)
axs[3].set_ylabel("PSD (Hz)", color="blue", fontsize=10)
axs[3].set_yscale('log')
axs[3].set_xlim([50, 600])
axs[3].set_ylim([10, 10000])

plt.grid()
plt.show()

"""
Polished
"""

# pgSFQQ1 = []
# l_i = 1
# for i in range(len(pgSFQ)):
#     pgSFQQ1.append(pgSFQ[i]*eQ1[i])
#
# pgSFQQ1_scaled = []
# p2QP = 8e-3 # photon to QP conversion rate
# base = 1000
# for i in range(len(pgSFQ)):
#     pgSFQQ1_scaled.append(pgSFQQ1[i]*eQ1[i]*p2QP+base)
#
# fig, axs = plt.subplots(2)
# axs_02 = axs[0].twinx()
#
# axs_02.plot(f_SFQ[l_i:], eSFQ[l_i:], color="red", marker="o", label='Radiator Efficiency')
# axs_02.plot(f_Q1[l_i:], eQ1[l_i:], color="yellow", marker="o", label='Receiver Efficiency')
# axs_02.set_ylabel('Coupling Efficiency', color="black", fontsize=10)
# axs_02.set_yscale('log')
# axs_02.set_ylim([1e-4, 1e-1])
# axs_02.legend(loc=1)
#
# axs[0].plot(f_Q1[l_i:], pgSFQQ1[l_i:], color="orange", marker="o", label='Photon Rate')
# axs[0].set_ylabel("Photons/Sec", color="black", fontsize=10)
# axs[0].set_xlim([50, 700])
# axs[0].set_ylim([1e5, 1e9])
# axs[0].set_yscale('log')
# axs[0].grid(True, which="both")
# axs[0].legend(loc=4)
# axs[0].set_xlabel("Antenna Frequency (GHz)", color="black", fontsize=10)
#
# axs_12 = axs[1].twinx()
# axs_12.plot(f_Q1[l_i:], pgSFQQ1_scaled[l_i:], color="orange", linestyle='--', label='Photon Rate Scaled')
# axs_12.set_ylabel("Scaled Photons Rate (Hz)", color="black", fontsize=10)
# axs_12.set_xlim([50, 700])
# axs_12.set_ylim([500, 10000])
# axs_12.set_yscale('log')
# axs_12.legend(loc=1)
#
# axs[1].plot(Q1_PSD[:, 0]*f_SIM, Q1_PSD[:, 1], color='k', label='Q1 Measurement 200fF')
# axs[1].set_xlabel("Radiator Josephson Frequency (GHz)", color="black", fontsize=10)
# axs[1].set_ylabel("Parity Rate (Hz)", color="black", fontsize=10)
# axs[1].set_yscale('log')
# axs[1].set_xlim([50, 700])
# axs[1].set_ylim([500, 10000])
# axs[1].grid(True, which="both")
# axs[1].legend(loc=4)
#
# plt.show()