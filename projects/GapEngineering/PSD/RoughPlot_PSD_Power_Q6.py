import matplotlib.pyplot as plt

power = [-20, -16, -13, -10]
f_Q6_Clean = [1000/0.4855, 1000/0.4705, 1000/0.47385, 1000/0.49975]
f_Q6_Dirty = [1000/0.46729, 1000/0.45107, 1000/0.44786,  1000/0.45201]
f_Q6_Base_Clean = [1000/0.49529, 1000/0.49529, 1000/0.49529, 1000/0.49529]
f_Q6_Base_Dirty = [1000/0.44417, 1000/0.44417, 1000/0.44417, 1000/0.44417]
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q6_Clean, '-+', label='Q6 Poison Clean')
plt.plot(power, f_Q6_Dirty, '-o', label='Q6 Poison Dirty')
plt.plot(power, f_Q6_Base_Clean, '--v', label='Q6 No Poison Clean')
plt.plot(power, f_Q6_Base_Dirty, '--s', label='Q6 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('PSD Knee Freq (Hz)')
plt.title('QP Tunneling Rate vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()