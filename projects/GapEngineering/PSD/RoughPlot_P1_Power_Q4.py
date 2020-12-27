import matplotlib.pyplot as plt

power = [-20, -16, -13, -10]
f_Q4_Clean = [0.03394, 0.03458, 0.03859, 0.03831]
f_Q4_Dirty = [0.04619, 0.04918, 0.07112,  0.1251]
f_Q4_Base_Clean = [0.02896, 0.02896, 0.02896, 0.02896]
f_Q4_Base_Dirty = [0.04385, 0.04385, 0.04385, 0.04385]
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q4_Clean, '-+', label='Q4 Poison Clean')
plt.plot(power, f_Q4_Dirty, '-o', label='Q4 Poison Dirty')
plt.plot(power, f_Q4_Base_Clean, '--v', label='Q4 No Poison Clean')
plt.plot(power, f_Q4_Base_Dirty, '--s', label='Q4 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('P1 Occupation')
plt.title('P1 Occupation vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()