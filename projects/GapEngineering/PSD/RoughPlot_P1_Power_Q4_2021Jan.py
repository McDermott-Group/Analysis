import matplotlib.pyplot as plt
import numpy as np

power = np.arange(-20, -5, 1)
P1_Q4_Clean = [0.032, 0.0306, 0.0322, 0.0305, 0.03088, 0.03048, 0.03118, 0.034794,
               0.03835, 0.03747, 0.0402, 0.0418, 0.04751, 0.05063, 0.0542]
P1_Q4_Dirty = [0.04574, 0.0446, 0.04293, 0.04583, 0.04412, 0.048852, 0.0502, 0.06504,
               0.06524, 0.07728, 0.0945, 0.128076, 0.1395, 0.16393, 0.1966]
# P1_Q4_Base_Clean = [0.0274]
P1_Q4_Base_Clean = np.ones(15) * 0.0274
# P1_Q4_Base_Dirty = [0.03267]
P1_Q4_Base_Dirty = np.ones(15) * 0.03267
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, P1_Q4_Clean, '-+', label='Q4 Poison Clean')
plt.plot(power, P1_Q4_Dirty, '-o', label='Q4 Poison Dirty')
plt.plot(power, P1_Q4_Base_Clean, '--v', label='Q4 No Poison Clean')
plt.plot(power, P1_Q4_Base_Dirty, '--s', label='Q4 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('P1 Occupation')
plt.title('P1 Occupation vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()