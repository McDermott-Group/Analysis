import matplotlib.pyplot as plt
import numpy as np

# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 =[0.014641618690277242, 0.11298288892161439, 0.0030341285261799786, 0.002874418539831116]
Neg20 =[0.013928752984883807, 0.11036252634721502, 0.003047926199988066, 0.00285447019016742]
Neg10 = [0.021226507135723483, 0.09725433509546787, 0.0030020584727337637, 0.0028129654165653137]
Neg6 =[0.02898589639299382, 0.15193047172811872, 0.0029352238867866847, 0.002579433885931115]
Neg5 = [0.0415029683460855, 0.11951930665421667, 0.0015669035073672497, 0.0020137493335408755]
Neg4 = [0.039375451168775506, 0.13486290571631557, 0.0025980536148426376, 0.0024815201751751366]
Neg3 = [0.034457829612526907, 0.1405104493578824, 0.0027480702029042234, 0.002509981370317287]
Neg2 = [0.034679277381682305, 0.1261449019071591, 0.0027909814431390087, 0.002457101503778163]
Neg1 = [0.03568966704829278, 0.13029382873812162, 0.0017767566677624807, 0.0017713759872487556]
Neg0 = [0.046936126500354926, 0.12181211764107111, 0.0016282017023108806, 0.0015980843389965034]

# Initialize data list
power_list = [Neg20, Neg10, Neg6, Neg5, Neg4, Neg3, Neg2, Neg1, Neg0]
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

power = [-20, -10, -6, -5, -4, -3, -2, -1, 0]


fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Relative Poison Power (dBm)')
ax1.set_ylabel('Parity Lifetime (ms)', color=color)
# ax1.set_ylim(1.5, 3.5)
ax1.plot(power, PSD_clean_list, '-+', label='Q6 Clean PSD', color=color)
ax1.plot(power, PSD_clean_base, '-', label='Q6 Clean PSD Base', color=color)
ax1.plot(power, PSD_dirty_list, '-o', label='Q6 Dirty PSD', color=color)
ax1.plot(power, PSD_dirty_base, '-', label='Q6 Dirty PSD Base', color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color ='tab:blue'
ax2.set_ylabel('P1', color=color)
# ax2.set_ylim(0, 1)
ax2.plot(power, P1_clean_list, '--+', label='Q6 Clean P1', color=color)
ax2.plot(power, P1_clean_base, '-', label='Q6 Clean P1 Base', color=color)
ax2.plot(power, P1_dirty_list, '--o', label='Q6 Dirty P1', color=color)
ax2.plot(power, P1_dirty_base, '-', label='Q6 Dirty P1 Base', color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.legend(loc=2)
plt.grid()
plt.show()

