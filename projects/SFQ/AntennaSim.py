from SFQlib import add_CST_SimResultFromTxt
import numpy as np
import matplotlib.pyplot as plt

# from scipy import fft
from numpy.fft import rfft, rfftfreq

i1_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsCurrentPort1.txt"
i2_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsCurrentPort1.txt"
# v1_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsVoltagePort1.txt"
# v2_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsVoltagePort2.txt"

v1_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsVoltagePort1_long.txt"
v2_file = "Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/CST_Sim/TimeDomainSignalsVoltagePort2_long.txt"

R1 = 1.0    # patch antenna resistance
R2 = 8.0e3  # circmon antenna resistance

t = add_CST_SimResultFromTxt(v1_file)[0]
# i1 = add_CST_SimResultFromTxt(i1_file)[1]
# i2 = add_CST_SimResultFromTxt(i2_file)[1]
v1 = add_CST_SimResultFromTxt(v1_file)[1]
v2 = add_CST_SimResultFromTxt(v2_file)[1]
print(len(v1))

p1 = 0.0
p2 = 0.0
p2_1 = 0.0
for i in range(len(v1)):
    if i <= 200:
        p1 = p1 + v1[i]**2/R1
    # if i <= 500:
    p2 = p2 + v2[i]**2/R2
    # if i > 500:


# print('p1=', p1)
# print('p2=', p2)
print('p2/p1=', p2/p1)
# print('p2_1/p1=', p2_1/p1)

# plt.plot(t[:400], v1[:400])
plt.plot(t, v1)
plt.plot(t, v2)
plt.show()

"""FFT"""
# n = len(t)
# yf = rfft(v2)
# data_step = (t[-1]-t[0])/(n-1)
# xf = rfftfreq(n, data_step)
# plt.plot(xf, np.abs(yf))
# plt.xlabel('THz')
# plt.ylabel('Amp(arb)')
# plt.show()

"""FFT Test"""
# t_test = np.linspace(0, 10, 10001)
# y_test = np.sin(5*2*np.pi*t_test)
# # plt.plot(t_test, y_test)
# yf = rfft(y_test)
# data_step = (t_test[-1]-t_test[0])/(10001.0)
# xf = rfftfreq(10001, data_step)
# plt.plot(xf, np.abs(yf))
# plt.show()
