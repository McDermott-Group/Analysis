import matplotlib.pyplot as plt
import numpy as np


# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]

# Initialize data list
power_list = [Neg20, Neg19, Neg18, Neg17, Neg16, Neg15, Neg14, Neg13, Neg12, Neg11, Neg10, Neg9, Neg8, Neg7, Neg6]
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

power = [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6]

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

ax2 = ax1.twinx()
color ='tab:blue'
ax2.set_ylabel('P1', color=color)
# ax2.set_ylim(0, 1)
ax2.plot(power, P1_clean_list, '--+', label='Q4 Clean P1', color=color)
ax2.plot(power, P1_clean_base, '-', label='Q4 Clean P1 Base', color=color)
ax2.plot(power, P1_dirty_list, '--o', label='Q4 Dirty P1', color=color)
ax2.plot(power, P1_dirty_base, '-', label='Q4 Dirty P1 Base', color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.legend(loc=2)
plt.grid()
plt.show()
