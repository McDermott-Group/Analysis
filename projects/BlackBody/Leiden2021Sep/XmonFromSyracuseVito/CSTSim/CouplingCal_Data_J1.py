from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np

### parameters to be tuned
e_eff = 6 # limit (1, 6.5), the voltage can also be built in to have a larger range
C_eff = 75*1e-21   # Commonly used (50-100)
Jbias_offset = 1    # mDAC should be +-1 mDAC basically +-5 GHz
k = 1   # Coupling between radiator and receiver, this could be larger than one due to the
        # fact we can generate QPs locally at the recevier's test pad

k1 = 0.0   # coupling between on chip phonon mitigation
f_SIM = 0.97
### parameteres tuned done

# JJ7 = [4.2*1e3, None, 0, 1000*150]   #[R, L, C, A]
JJ1 = [33*1e3, None, 0, 180*150]   #[R, L, C, A]
JQ1 = [15.67*1e3, None, 0, 360*150]

fileJ1 = "Z_syracuse_injector_with_leads.txt"
fileQ1 = "Z_syracuse_xmon_5.1GHz_no_Coupler.txt"

J1 = AntennaCoupling()
J1.import_data(fileJ1, JJ1)
f = J1.Antenna["f"]
ecJ1 = J1.e_c_dB
pgJ1 = J1.p_g
refJ1 = J1.ref

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1)
f_Q1 = Q1.Antenna["f"]
ecQ1 = Q1.e_c_dB
eQ1 = Q1.e_c
Z_ReQ1 = Q1.Antenna["Z_Re"]

Q3_PSD = np.array([
    [0, 106.35], [10, 104.29], [20, 104.98], [30, 108.07], [40, 107.02],
    [50, 104.26], [60, 103.02], [70, 100.7], [80, 106.52], [90, 106.48],
    [100, 99.7], [105, 102.78], [110, 113.84], [115, 102.97], [120, 98.58],
    [125, 108.68], [130, 108.93], [135, 110.95], [140, 112.91], [145, 106.57],
    [150, 110.2], [155, 112.03], [160, 110.75], [165, 115.04], [170, 111.42],
    [175, 129.16], [180, 129.86], [185, 134.93], [190, 129.63], [195, 208.11],
    [200, 168.99], [205, 155.97], [210, 151.51], [215, 127.42], [220, 140.36],
    [225, 133.45], [230, 128.69], [235, 139.76], [240, 138.28], [245, 155.1],
    [250, 169.76], [255, 179.9], [260, 200.54], [265, 340.67], [270, 377.06],
    [275, 558.06], [280, 570.87], [285, 549.9], [290, 483.57], [295, 355.19],
    [300, 219.61], [305, 254.13], [310, 209.86], [315, 197.07], [320, 165.21],
    [325, 157.97], [330, 151.87], [335, 148.3], [340, 158.34], [345, 131.75],
    [350, 135.81], [355, 140.82], [360, 133.66], [365, 132.2], [370, 128.49],
    [375, 133.89], [380, 136.2], [385, 146.84], [390, 132.96], [395, 124.15],
    [400, 131.54], [405, 126.31], [410, 128.52], [415, 131.69], [420, 144.64],
    [425, 137.99], [430, 144.86], [435, 153.82], [440, 143.3], [445, 139.44],
    [450, 144.97], [455, 146.19], [460, 141.58], [465, 148.59], [470, 143.99],
    [475, 140.83], [480, 164.62], [485, 145.94], [490, 166.24], [495, 163.89],
    [500, 168.67], [505, 168.55], [510, 158.79], [515, 156.61], [520, 162.74],
    [525, 158.32], [530, 164.04], [535, 160.66], [540, 163.78], [545, 171.71],
    [550, 166.11], [555, 177.25], [560, 176.65], [565, 185.42], [570, 195.07],
    [575, 194.23], [580, 193.45], [585, 199.59], [590, 215.71], [595, 218.82],
    [600, 215.64], [605, 236.99], [610, 240.65], [615, 240.53], [620, 247.78],
    [625, 254.46], [630, 256.9], [635, 269.65], [640, 274.44], [645, 279.23],
    [650, 297.19], [655, 289.41], [660, 303.32], [665, 305.88], [670, 313.94],
    [675, 319.42], [680, 331.56], [685, 336.2], [690, 346.79], [695, 350.34]
])

Q3_PSD[:, 0] = (Q3_PSD[:, 0])

f_scale = np.sqrt(e_eff/6.0)
f = f/1e9
f = f*f_scale
f_Q1 = f_Q1/1e9
f_Q1 = f_Q1*f_scale

### dB power
# fig, axs = plt.subplots(4)
# axs[0].plot(f, ecJ1, color="black", marker="o")
# axs[0].set_ylabel("Radiator (dB)", color="black", fontsize=10)
# axs[0].set_xlim([50, 600])
# axs[0].set_ylim([-30, -10])

# axs[1].plot(f_Q1, ecQ1, color="green", marker="o")
# axs[1].set_ylabel("Receiver QB (dB)", color="green", fontsize=10)
# axs[1].set_xlim([50, 600])
# axs[1].set_ylim([-50, 0])
#
# axs[2].plot(f_Q1, ecQ1+ecJ1*k, color="red", marker="o")
# axs[2].set_ylabel("Sum (dB)", color="red", fontsize=10)
# axs[2].set_xlim([50, 600])
# axs[2].set_ylim([-70, -20])
#
# axs[3].plot(Q3_PSD[:, 0]*f_SIM, Q3_PSD[:, 1], color='k', label='Q3_PSD')
# axs[3].set_xlabel("Freq (GHz)", color="black", fontsize=10)
# axs[3].set_ylabel("PSD (Hz)", color="blue", fontsize=10)
# axs[3].set_yscale('log')
# axs[3].set_xlim([50, 600])
# axs[3].set_ylim([100, 1000])
#
# plt.grid()
# plt.show()

### Absolute photon rate
fig, axs = plt.subplots(4)
axs[0].plot(f, pgJ1, color="black", marker="o")
axs[0].set_ylabel("Photons/Sec", color="black", fontsize=10)
axs[0].set_xlim([50, 600])
axs[0].set_ylim([1e6, 1e10])
axs[0].set_yscale('log')
axs[0].grid(True, which="both")

axs[1].plot(f_Q1[100:], eQ1[100:], color="green", marker="o")
axs[1].set_ylabel("Receiver QB Gamma", color="green", fontsize=10)
axs[1].set_xlim([50, 600])
axs[1].set_ylim([1e-5, 1e0])
axs[1].set_yscale('log')
axs[1].grid(True, which="both")


pgJ1Q = []
for i in range(len(pgJ1)):
    pgJ1Q.append(pgJ1[i]*eQ1[i]*refJ1[i])
# print(pgJ1Q)
axs[2].plot(f_Q1[100:], pgJ1Q[100:], color="red", marker="o")
# axs[2].plot(f_Q1[100:], refJ1[100:], color="red", marker="o")
axs[2].set_ylabel("Photons/Sec", color="red", fontsize=10)
axs[2].set_xlim([50, 600])
# axs[2].set_ylim([1e3, 1e8])
axs[2].set_yscale('log')
axs[2].grid(True, which="both")

axs[3].plot(Q3_PSD[:, 0]*f_SIM, Q3_PSD[:, 1], color='k', label='Q3_PSD')
axs[3].set_xlabel("Freq (GHz)", color="black", fontsize=10)
axs[3].set_ylabel("PSD (Hz)", color="blue", fontsize=10)
axs[3].set_yscale('log')
axs[3].set_xlim([50, 600])
axs[3].set_ylim([100, 1000])
axs[3].grid(True, which="both")


plt.show()