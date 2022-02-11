from antennalib import QP_Up, Up_array
import numpy as np
import matplotlib.pyplot as plt

#Bias = 32.2
GammaUp=624.5424275627762
time=np.array([ 4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15.,])
occ=np.array([0.01244255, 0.01309302, 0.01369066, 0.01434748, 0.01497818, 0.01581184, 0.01708013, 0.01663164, 0.017047, 0.01783819, 0.01867434, 0.01968149])
error=np.array([0.00046377, 0.00037788, 0.00059402, 0.00048323, 0.00068067, 0.00065523, 0.00061548, 0.0005248, 0.00048127, 0.0005305, 0.00063201, 0.00062333])
occ_fit=np.array([0.01250806, 0.0131326, 0.01375715, 0.01438169, 0.01500623, 0.01563077, 0.01625531, 0.01687986, 0.0175044, 0.01812894, 0.01875348, 0.01937803])

plt.plot(list([0])+list(time),list([-(time[0])*(10**-6)*GammaUp+occ_fit[0]])+list(occ_fit))
plt.errorbar(time, occ, yerr=error,fmt='o', label='Q4 GammaUp={0:.4g} Hz'.format(GammaUp))
plt.xlim([0,15.5])
plt.xlabel('time (us)')
plt.ylabel('P1')
plt.grid()
plt.legend()
plt.show()

#Bias = 28.3
GammaUp=543.6357980508438
time=np.array([ 4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15.,])
occ=np.array([0.01157711, 0.01223969, 0.0131538, 0.0132512, 0.01401482, 0.01507381, 0.01504554, 0.01587637, 0.01558827, 0.01710131, 0.01662816, 0.01804158])
error=np.array([0.00060196, 0.00045671, 0.00050455, 0.00065888, 0.00057993, 0.00063427, 0.00068612, 0.00060045, 0.00084991, 0.0006179, 0.00061,0.00068525])
occ_fit=np.array([0.01180931, 0.01235294, 0.01289658, 0.01344021, 0.01398385, 0.01452749, 0.01507112, 0.01561476, 0.01615839, 0.01670203, 0.01724567, 0.0177893 ])

plt.plot(list([0])+list(time),list([-(time[0])*(10**-6)*GammaUp+occ_fit[0]])+list(occ_fit))
plt.errorbar(time, occ, yerr=error,fmt='o', label='Q4 GammaUp={0:.4g} Hz'.format(GammaUp))
plt.xlim([0,15.5])
plt.xlabel('time (us)')
plt.ylabel('P1')
plt.grid()
plt.legend()
plt.show()
