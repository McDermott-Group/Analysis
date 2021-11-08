"""
PSD for Q1, Q2, Q3 at different J1 Bias with Sim slow bias
Data:
Z:\mcdermott-group\data\Antenna\Circmon\Liu\CW20180514A_Ox2\2021Nov02_QB124_PSD_J1Radiator_Sim
Fitting method Chris' no white noise verison

"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np




# Q3 = np.array([
#     [-280, 507.66], [-270, 552.25], [-260, 192.2], [-250, 200.82], [-240, 169.66],
# [240, 135.77], [250, 147.43], [260, 219.18], [270, 213.15], [280, 402.44]
# ])

Q3 = np.array([
    [240, 131.51], [250, 168.97], [260, 209.19], [270, 382.35], [280, 592.34],
    [290, 482.65], [300, 229.36]
])


### J1 Bias Conversion
Q3[:, 0] = (Q3[:, 0])

# plt.plot(Q3[:, 0], Q3[:, 1], color='y', label='Q3')
# plt.xlabel('Radiator Josephson Frequency (mDAC)')
# plt.ylabel('PSD (Hz)')
# plt.yscale('log')
# plt.grid()
# plt.legend(loc=1)
# plt.show()

f = 1
# f = 0.97
Al_gap = 380e-6
DAC_Al = 1e5*Al_gap/0.200
# plt.plot(Q1[::2, 0]*f, Q1[::2, 1], color='b', label='Q1')
# plt.plot(Q2[::2, 0]*f, Q2[::2, 1], color='r', label='Q2')
# plt.plot(Q3[::2, 0]*f, Q3[::2, 1], color='y', label='Q3')
plt.plot(Q3[:, 0]*f, Q3[:, 1], color='y', label='Q3')
# plt.plot(Q3[:, 0]*f, -Q3[:, 1], color='y', label='Q3')

# plt.plot(Q1[:, 0]*f, Q1[:, 1]*0.01, '-', label='Q2FromQ1')
# plt.plot(Q3[:, 0]*f, Q3[:, 1]*0.04, '-', label='Q2FromQ3')
# plt.plot(Q2[:, 0]*f, Q3[:, 1]*0.04+Q1[:, 1]*0.008, '--', label='Q2FromQ14')
# plt.plot(Q2[:, 0]*f, Q2[:, 1]-(Q3[:, 1]*0.04+Q1[:, 1]*0.008), '--', label='Q2-Q14')

plt.axvline(x=DAC_Al * f, color='k', linestyle='--', linewidth=4, label='JJ Al Gap')

plt.xlabel('Radiator Josephson Frequency (GHz)')
plt.ylabel('PSD (Hz)')
plt.yscale('log')
plt.grid()
plt.legend(loc=1)
# plt.xlim([0, 1500])
# plt.ylim([10, 100000])
plt.show()



