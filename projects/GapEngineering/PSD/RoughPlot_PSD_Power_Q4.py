import matplotlib.pyplot as plt

power = [-20, -16, -13, -10]
f_Q4_Clean = [342, 324, 322, 326]
f_Q4_Dirty = [361, 344, 380, 418]
f_Q4_Base_Clean = [320, 320, 320, 320]
f_Q4_Base_Dirty = [300, 300, 300, 300]
# f_Q6 = [242, 233, 235, 235]
# f_Q6_Base = [256, 256, 256, 256]
# '+', 'o', 'v', 's'
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q4_Clean, '-+', label='Q4 Poison Clean')
plt.plot(power, f_Q4_Dirty, '-o', label='Q4 Poison Dirty')
plt.plot(power, f_Q4_Base_Clean, '--v', label='Q4 No Poison Clean')
plt.plot(power, f_Q4_Base_Dirty, '--s', label='Q4 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('PSD Knee Freq (Hz)')
plt.title('QP Tunneling Rate vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()