"""
PSD for Q12 at different J7 Bias
Data:
Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\PSD2021Oct25_Q2WithJ7
Fitting method Chris' no white noise verison

"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np

Q1 = np.array([
    [0.0, 133.64], [0.005, 142.46], [0.01, 130.96], [0.015, 44439.41], [0.02, 629.71],
    [0.025, 1926.76], [0.0275, 8044.96], [0.03, 6184.95]
])

Q2 = np.array([
    [0.0, 155.45], [0.005, 157.37], [0.01, 155.79], [0.015, 165.68], [0.02, 445.68],
    [0.025, 1468.44], [0.0275, 1355.17], [0.03, 6121.79], [0.04, 6477.23]
])

"""
Q2 no good fit
 [0.045, 74757.51], [0.0499, 76559.81], [0.0549, 80510.31], [0.0599, 82672.48],
    [0.0649, 82651.46], [0.0699, 79224.2], [0.0749, 70857.6], [0.0799, 78456.67],
    [0.0849, 75701.37], [0.0899, 80643.0], [0.0949, 78892.78], [0.0999, 76166.65]
"""

### J7 Bias Conversion
Q1[:, 0] = (Q1[:, 0])*1000
Q2[:, 0] = (Q2[:, 0])*1000


# plt.plot(Q1[:, 0], Q1[:, 1], color='b', label='Q1')
# plt.plot(Q2[:, 0], Q2[:, 1], color='r', label='Q2')
# plt.plot(Q3[:, 0], Q3[:, 1], color='y', label='Q3')
# plt.xlabel('Radiator Josephson Frequency (mDAC)')
# plt.ylabel('PSD (Hz)')
# plt.yscale('log')
# plt.grid()
# plt.legend(loc=1)
# plt.show()

f = 4.604
# f = 1
Al_gap = 380e-6
DAC_Al = 1e5*Al_gap/0.952
# plt.plot(Q1[::2, 0]*f, Q1[::2, 1], color='b', label='Q1')
# plt.plot(Q2[::2, 0]*f, Q2[::2, 1], color='r', label='Q2')
# plt.plot(Q3[::2, 0]*f, Q3[::2, 1], color='y', label='Q3')
# plt.plot(Q1[:, 0]*f, Q1[:, 1], color='b', label='Q1')
plt.plot(Q2[:, 0]*f, Q2[:, 1], color='r', label='Q2')
# plt.plot(Q3[:, 0]*f, Q3[:, 1], color='y', label='Q3')

plt.axvline(x=DAC_Al * f, color='k', linestyle='--', linewidth=4, label='JJ Al Gap')

plt.xlabel('Radiator Josephson Frequency (GHz) Same Chip')
plt.ylabel('PSD (Hz)')
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.legend(loc=1)
plt.xlim([0, 800])
# plt.ylim([50, 20000])
plt.show()



