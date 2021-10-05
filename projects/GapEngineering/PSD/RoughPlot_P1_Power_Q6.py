import matplotlib.pyplot as plt

power = [-20, -16, -13, -10]
f_Q6_Clean = [0.0632, 0.0642, 0.0678, 0.0704]
f_Q6_Dirty = [0.1221, 0.1251, 0.1459, 0.1248]
f_Q6_Base_Clean = [0.05914, 0.05914, 0.05914, 0.05914]
f_Q6_Base_Dirty = [0.08576, 0.08576, 0.08576, 0.08576]
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q6_Clean, '-+', label='Q6 Poison Clean')
plt.plot(power, f_Q6_Dirty, '-o', label='Q6 Poison Dirty')
plt.plot(power, f_Q6_Base_Clean, '--v', label='Q6 No Poison Clean')
plt.plot(power, f_Q6_Base_Dirty, '--s', label='Q6 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('P1 Occupation')
plt.title('P1 Occupation vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()