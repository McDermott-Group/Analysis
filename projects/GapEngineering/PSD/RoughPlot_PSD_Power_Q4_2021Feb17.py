import matplotlib.pyplot as plt
import numpy as np


# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 = [0.018530865099552565, 0.07201861242659101, 0.0019956553520469556, 0.0018750782951790271]
Neg20 = [0.020585912308760876, 0.08705324513922055, 0.0020169592643044343, 0.0018040558920183031]
Neg19 = [0.020485699836699008, 0.07559138342078968, 0.0020051690522098855, 0.0018582262068804647]
Neg18 = [0.021496346143509152, 0.12232174369593937, 0.002005731640621227, 0.0018115412399837348]
Neg17 = [0.021112111253508616, 0.12077237388415502, 0.0020169474411580054, 0.0018422819300343387]
Neg16 = [0.02180844852409245, 0.1010919558884891, 0.0020208663956347806, 0.0018411593671858297]
Neg15 = [0.02209441001328936, 0.1312178618868927, 0.0020403593602501624, 0.0016976426612857255]
Neg14 = [0.021509901131529893, 0.09608901027393081, 0.0020811758251168596, 0.001888529147212627]
Neg13 = [0.022902869440556493, 0.08711945363077414, 0.002060626243945416, 0.001881427798010267]
Neg12 = [0.02356703488318083, 0.11539325601901353, 0.0020137739334743516, 0.0018240762654737636]
Neg11 = [0.02460260904459349, 0.09231133175172049, 0.002030483055091358, 0.0017841879013880825]
Neg10 = [0.025430769090984087, 0.08653970025305537, 0.0019945587512915615, 0.001769721308436739]
Neg9 = [0.026963456331042718, 0.07416659146140495, 0.0019375323969103237, 0.001819208856116842]
Neg8 = [0.027991057692758278, 0.09933187641514588, 0.001971204056054535, 0.0017817869470029251]
Neg7 = [0.03332933356140834, 0.07362910679957971, 0.0017688928591725278, 0.0017126931019757247]
Neg6 = [0.034354041297449164, 0.0736866585445068, 0.0016677623315069323, 0.0016270699691382747]

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
# # color = 'tab:red'
# color = 'tab:red'
# ax1.set_xlabel('Relative Poison Power (dBm)')
# ax1.set_ylabel('Parity Lifetime (ms)', color=color)
# # ax1.set_ylim(1.5, 3.5)
# ax1.plot(power, PSD_clean_list, '-+', label='Q4 Clean PSD', color=color)
# ax1.plot(power, PSD_clean_base, '-', label='Q4 Clean PSD Base', color=color)
# ax1.plot(power, PSD_dirty_list, '-o', label='Q4 Dirty PSD', color=color)
# ax1.plot(power, PSD_dirty_base, '-', label='Q4 Dirty PSD Base', color=color)
# ax1.tick_params(axis='y', labelcolor=color)

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
plt.title('Parity Lifetime vs Poison Power')
plt.legend(loc=1)
plt.rc('font', size=22)
plt.grid()
plt.show()
