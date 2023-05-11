import matplotlib.pyplot as plt
import numpy as np

label_font = 18
tick_font = 16
legend_font = 14

day = np.array([1, 2, 3, 4, 14, 15])
day_fit = np.arange(1, 15.1, 0.1)
Q1 = np.array([23.2, 18.9, 18.2, 14.6, 10, 10])
Q2 = np.array([30.8, 27.6, 23.5, 18.8, 11.7, 12.0])
Q3 = np.array([27.3, 22.6, 20.2, 16.8, 10.3, 11.1])

Q1_params = np.polyfit(day, np.log(Q1), 1)
Q1_fit = np.exp(Q1_params[1])*np.exp(day_fit*Q1_params[0])

Q2_params = np.polyfit(day, np.log(Q2), 1)
Q2_fit = np.exp(Q2_params[1])*np.exp(day_fit*Q2_params[0])

Q3_params = np.polyfit(day, np.log(Q3), 1)
Q3_fit = np.exp(Q3_params[1])*np.exp(day_fit*Q3_params[0])

print(Q1_params[0], Q2_params[0], Q3_params[0])

plt.figure(figsize=(9, 6))
plt.plot(day_fit, Q1_fit, c='r')
plt.scatter(day, Q1, color='r')
plt.plot(day_fit, Q2_fit, c='k')
plt.scatter(day, Q2, color='k')
plt.plot(day_fit, Q3_fit, c='b')
plt.scatter(day, Q3, color='b')
plt.xlabel('Days after condensation', fontsize=label_font)
plt.ylabel('Charge parity switch rate ($s^{-1}$)', fontsize=label_font)
plt.ylim([5, 30])
plt.legend(frameon=False, loc=2, prop={'size': 14})
plt.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)

path = 'Z:\mcdermott-group\users\ChuanghongVincentLiu\Thesis\plotsNotInAntennaSFQPapers'
plt.savefig(path + '\QPThermalisation.pdf', format='pdf', bbox_inches='tight', dpi=900)

plt.show()