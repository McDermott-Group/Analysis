import matplotlib.pyplot as plt
import numpy as np

# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 = [0.01088, 0.11223, 0.0032682, 0.0030714]
Neg20 = [0.0117, 0.163234, 0.003305995, 0.0030117]
Neg19 = [0.01205445, 0.13125, 0.0033016, 0.0031228]
Neg18 = [0.0122843, 0.153764665, 0.00336, 0.0030478]
Neg17 = [0.0126, 0.1511, 0.00338, 0.00295]
Neg16 = [0.014697, 0.17147, 0.00334986, 0.00313996]
Neg15 = [0.0135014, 0.20299277, 0.00337222, 0.0030034]
Neg14 = [0.01409983, 0.135112, 0.0033582, 0.0030353]
Neg13 = [0.014199, 0.1678699, 0.0033385, 0.002923369]
Neg12 = [0.017434, 0.13180, 0.003337, 0.003108686]
Neg11 = [0.016272, 0.20571894, 0.0032894, 0.0029882]
Neg10 = [0.01881, 0.16950, 0.0031829, 0.0029656]
Neg9 = [0.02114436, 0.1204981, 0.00309277, 0.0029963]
Neg8 = [0.0227389, 0.1171387, 0.00316297, 0.0029226]
Neg7 = [0.0234465, 0.1220778, 0.00310084, 0.002934545]
Neg6 = [0.023717, 0.1345021, 0.0031196697, 0.00289809]

# Initialize data list
power_list = [Neg20, Neg19, Neg18, Neg17, Neg16, Neg15, Neg14, Neg13, Neg12, Neg11, Neg10, Neg9, Neg8, Neg7, Neg6]
P1_clean_list = np.zeros(len(power_list))
P1_dirty_list = np.zeros(len(power_list))
PSD_clean_list = np.zeros(len(power_list))
PSD_dirty_list = np.zeros(len(power_list))

P1_clean_base = np.ones(len(power_list)) * Neg100[0]
P1_dirty_base = np.ones(len(power_list)) * Neg100[1]
PSD_clean_base = np.ones(len(power_list)) * Neg100[2]*1000
PSD_dirty_base = np.ones(len(power_list)) * Neg100[3]*1000

# Update data
for i, power in enumerate(power_list):
    P1_clean_list[i]=power[0]
    P1_dirty_list[i]=power[1]
    PSD_clean_list[i]=power[2]*1000
    PSD_dirty_list[i]=power[3]*1000

power = [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6]


fig, ax1 = plt.subplots()
# color = 'tab:red'
color = 'tab:blue'
ax1.set_xlabel('Relative Poison Power (dBm)')
ax1.set_ylabel('Parity Lifetime (ms)', color=color)
# ax1.set_ylim(1.5, 3.5)
ax1.plot(power, PSD_clean_list, '-+', label='Q6 Clean PSD', color=color)
ax1.plot(power, PSD_clean_base, '-', label='Q6 Clean PSD Base', color=color)
ax1.plot(power, PSD_dirty_list, '-o', label='Q6 Dirty PSD', color=color)
ax1.plot(power, PSD_dirty_base, '-', label='Q6 Dirty PSD Base', color=color)
ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()
# color ='tab:blue'
# ax2.set_ylabel('P1', color=color)
# # ax2.set_ylim(0, 1)
# ax2.plot(power, P1_clean_list, '--+', label='Q6 Clean P1', color=color)
# ax2.plot(power, P1_clean_base, '-', label='Q6 Clean P1 Base', color=color)
# ax2.plot(power, P1_dirty_list, '--o', label='Q6 Dirty P1', color=color)
# ax2.plot(power, P1_dirty_base, '-', label='Q6 Dirty P1 Base', color=color)
# ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.title('Parity Lifetime vs Poison Power')
plt.legend(loc=1)
plt.rc('font', size=22)
plt.grid()
plt.show()

