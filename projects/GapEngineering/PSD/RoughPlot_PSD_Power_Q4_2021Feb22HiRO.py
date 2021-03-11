import matplotlib.pyplot as plt
import numpy as np


# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 = [0.030541418092291963, 0.1427882530013499, 0.0015700730235426084, 0.001344505754424018]
Neg20 = [0.030743663956825203, 0.0835101089394041, 0.0015599272427418709, 0.0014182503823829302]
Neg17 = [0.031117656538502344, 0.08850083608755888, 0.0015262675159915657, 0.0014757368674501072]
Neg14 = [0.03613660434538867, 0.11072989159934954, 0.0015412394027753166, 0.0014391777835268643]
Neg11 = [0.05452062106278037, 0.07338136586494447, 0.0013990939820973133, 0.0014086569102117203]
Neg8 = [0.03826212291913469, 0.10769252170209713, 0.0015222648468499733, 0.0013747875499869207]

# Initialize data list
power_list = [Neg20, Neg17,Neg14,Neg11, Neg8]
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

power = [-20, -17, -14, -11, -8]

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
