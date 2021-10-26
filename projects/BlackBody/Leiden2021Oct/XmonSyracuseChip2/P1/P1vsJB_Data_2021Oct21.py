"""
P1 for Q1, Q2, Q3 at different J7 Bias
Data:


"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np


Q1 = np.array([])

### J7 Bias Conversion
Q1[:, 0] = (Q1[:, 0])*1000
Q2[:, 0] = (Q2[:, 0])*1000
Q4[:, 0] = (Q4[:, 0])*1000


f = 4.604
Al_gap = 380e-6
DAC_Al = 1e5*Al_gap/0.952
plt.errorbar(Q1[:, 0]*f, Q1[:, 1], yerr=Q1[:, 2]/np.sqrt(150), color='b', label='Q1')
plt.errorbar(Q2[:, 0]*f, Q2[:, 1], yerr=Q2[:, 2]/np.sqrt(150), color='r', label='Q2')
plt.errorbar(Q4[:, 0]*f, Q4[:, 1], yerr=Q4[:, 2]/np.sqrt(150), color='y', label='Q4')

plt.axvline(x=DAC_Al * f, color='k', linestyle='--', linewidth=4, label='JJ Al Gap')

plt.xlabel('Radiator Josephson Frequency (GHz)')
plt.ylabel('P1')
plt.yscale('log')
plt.grid()
plt.legend(loc=1)
# plt.xlim([0, 1500])
# plt.ylim([10, 100000])
plt.show()



