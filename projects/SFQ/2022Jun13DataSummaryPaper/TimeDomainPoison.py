# from SFQlib import DeepSubharmonics
import matplotlib.pyplot as plt
import numpy as np

# file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')
#
# date = '2022Jun28PoisonTimeDomainRealV1'
# # file_Number = np.arange(0, 300, 1)
# file_Number = np.arange(0, 150, 1)
# # experiment_name = ('PoisonQ1')
# experiment_name = ('PoisonQ2')
# file_list_Q1 = [file_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
# P1 = DeepSubharmonics()
# P1.add_data_from_matlab(file_list_Q1)
# P1.plot()
# #
# t = P1.Time_Dep
# P1_avg = P1.occ_1D_avg
# P1_se = P1.occ_1D_std
# P1_data = []
# for i in range(len(t)):
#     p1_data = [t[i], P1_avg[i], P1_se[i]]
#     P1_data.append(p1_data)
# np.savetxt('Q2P1.txt', P1_data)

# mpl.rc('font', family='Arial')
label_font = 18
tick_font = 16
legend_font = 14

Q1File = 'Q1P1.txt'
Q2File = 'Q2P1.txt'
Q1_t = np.loadtxt(Q1File, usecols=[0], skiprows=3)
Q1_P1_avg = np.loadtxt(Q1File, usecols=[1], skiprows=3)
Q1_P1_std = np.loadtxt(Q1File, usecols=[2], skiprows=3)
Q2_P1_avg = np.loadtxt(Q2File, usecols=[1], skiprows=3)
Q2_P1_std = np.loadtxt(Q2File, usecols=[2], skiprows=3)

plt.errorbar(Q1_t/1e3, Q1_P1_avg*100, yerr=Q1_P1_std*100, color='r', fmt="o", capsize=3.0, ms=5.0)
plt.errorbar(Q1_t/1e3, Q2_P1_avg*100, yerr=Q2_P1_std*100, color='b', fmt="o", capsize=3.0, ms=5.0)
plt.xlabel('Idle ($\mu$s)', fontsize=label_font)
plt.ylabel('$P_{1} (\%)$', fontsize=label_font)
plt.legend(frameon=False, loc=2, prop={'size': 14})
plt.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)

path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
plt.savefig(path + '\PoisonP1.pdf', format='pdf', bbox_inches='tight', dpi=900)

plt.show()
