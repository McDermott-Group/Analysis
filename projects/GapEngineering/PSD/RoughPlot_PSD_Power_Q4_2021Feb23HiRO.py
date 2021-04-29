import matplotlib.pyplot as plt
import numpy as np


# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 = [0.024157771567884293, 0.05060142337283608, 0.0015355092694927504, 0.0013987002118677342]
Neg20 = [0.02682339040906014, 0.06791523867532573, 0.0015316934279936782, 0.0015431016745079723]
Neg17 = [0.03222608331693518, 0.07198309784295098, 0.0015077936827345877, 0.0013811185943356026]
Neg14 = [0.05767647219634807, 0.08177826550662705, 0.001449342789731053, 0.001398161872708062]
Neg11 = [0.04762068201807486, 0.1289378015052424, 0.0015418329397289148, 0.0012924942981941796]

# Initialize data list
power_list = [Neg20, Neg17, Neg14, Neg11]
P1_clean_list = np.zeros(len(power_list))
P1_dirty_list = np.zeros(len(power_list))
PSD_clean_list = np.zeros(len(power_list))
PSD_dirty_list = np.zeros(len(power_list))

# Update data
for i, power in enumerate(power_list):
    P1_clean_list[i]=power[0]
    P1_dirty_list[i]=power[1]
    PSD_clean_list[i]=power[2]*1000
    PSD_dirty_list[i]=power[3]*1000

P1_clean_base = np.ones(len(power_list)) * Neg100[0]
P1_dirty_base = np.ones(len(power_list)) * Neg100[1]
PSD_clean_base = np.ones(len(power_list)) * Neg100[2]*1000
PSD_dirty_base = np.ones(len(power_list)) * Neg100[3]*1000

power = [-20, -17, -14, -11]

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Relative Poison Power (dBm)')
ax1.set_ylabel('Parity Lifetime (ms)', color=color)
# ax1.set_ylim(1.5, 3.5)
ax1.plot(power, PSD_clean_list, '-+', label='Q4 Clean PSD', color=color)
ax1.plot(power, PSD_clean_base, '-', label='Q4 Clean PSD Base', color=color)
ax1.plot(power, PSD_dirty_list, '-o', label='Q4 Dirty PSD', color=color)
ax1.plot(power, PSD_dirty_base, '-', label='Q4 Dirty PSD Base', color=color)
ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()
# color ='tab:blue'
# ax2.set_ylabel('P1', color=color)
# # ax2.set_ylim(0, 1)
# ax2.plot(power, P1_clean_list, '--+', label='Q4 Clean P1', color=color)
# # ax2.plot(power, P1_clean_base, '-', label='Q4 Clean P1 Base', color=color)
# ax2.plot(power, P1_dirty_list, '--o', label='Q4 Dirty P1', color=color)
# # ax2.plot(power, P1_dirty_base, '-', label='Q4 Dirty P1 Base', color=color)
# ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.legend(loc=2)
plt.grid()
plt.show()
