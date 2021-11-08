from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np
JSFQ = [16*1e3, None, 0, 100*200]   #[R, L, C, A]
JQ1 = [16.6*1e3, 18.3*1e-9, 0, 193.8*121.8]   #
JQ2 = [13.2*1e3, 14.6*1e-9, 0, 350*126.6]   #
JQ3 = [19*1e3, 21*1e-9, 0, 310*126]   #    some issue with Q3
JQ4 = [15*1e3, 19.9*1e-9, 0, 184.4*122.5*2]   #    some issue with Q3

fileSFQ = "SFQ_JJRadiator.txt"
# fileQ1 = "Q1.txt"
fileQ1 = "Q1_leads.txt"
fileQ2 = "Q2.txt"
fileQ3 = "Q3.txt"
fileQ4 = "Q4.txt"

SFQ = AntennaCoupling()
SFQ.import_data(fileSFQ, JSFQ)
f_SFQ = SFQ.Antenna["f"]
ecSFQ = SFQ.e_c_dB
# SFQ.plot()

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1)
# Q1.plot()
ecQ1 = Q1.e_c_dB
Z_ReQ1 = Q1.Antenna["Z_Re"]
f = Q1.Antenna["f"]

Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2)
ecQ2 = Q2.e_c_dB
Z_ReQ2 = Q2.Antenna["Z_Re"]

Q3 = AntennaCoupling()
Q3.import_data(fileQ3, JQ3)
ecQ3 = Q3.e_c_dB
Z_ReQ3 = Q3.Antenna["Z_Re"]

Q4 = AntennaCoupling()
Q4.import_data(fileQ4, JQ4)
ecQ4 = Q4.e_c_dB
Z_ReQ4 = Q4.Antenna["Z_Re"]

ecSFQ = ecSFQ*1.0
k = 1

# print('len(f)=', len(f))
# print('sum(f)=', np.sum(f))
# # print('f[0]=', f[0])
# # print('f[0:10]=', f[0:10])
# # print('f[-1]=', f[-1])
# # print('len(ecQ1)=', len(ecQ1))
# # print('len(ecQ2)=', len(ecQ2))
# # print('len(ecQ4)=', len(ecQ4))
# print('len(f_SFQ)=', len(f_SFQ))
# print('sum(f_SFQ)=', np.sum(f_SFQ))
# print('sum(f)-sum(f_SFQ)=', np.sum(f)-np.sum(f_SFQ))
# # print('f_SFQ[0]=', f_SFQ[0])
# # print('f_SFQ[0:10]=', f_SFQ[0:10])
# # print('f_SFQ[-1]=', f_SFQ[-1])
# # print('len(ecSFQ)=', len(ecSFQ))
# diff = np.setdiff1d(f, f_SFQ)
# print('diff=', diff)

plt.plot(f, ecQ1+ecSFQ*k, color='b', label='Q1')
plt.plot(f, ecQ2+ecSFQ*k, color='r', label='Q2')
# plt.plot(f, ecQ3+ecSFQ*k, color='k', label='Q3')
plt.plot(f, ecQ4+ecSFQ*k, color='y', label='Q4')
plt.plot(f_SFQ, ecSFQ, color='y', label='SFQ')
plt.xlabel('Freq')
plt.ylabel('Coupling efficiency')
plt.legend()
plt.ylim([-40, -10])
plt.show()

# plt.plot(f, Z_ReQ1, color='b', label='Q1')
# plt.plot(f, Z_ReQ2, color='r', label='Q2')
# plt.plot(f, Z_ReQ3, color='k', label='Q3')
# plt.plot(f, Z_ReQ4, color='y', label='Q4')
# # plt.plot(f, ecSFQ, color='y', label='SFQ')
# plt.xlabel('Freq')
# plt.ylabel('Coupling efficiency')
# plt.legend()
# plt.show()
