import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Q4Corr\General\09-10-20\T1\MATLABData\T1_{:03d}.mat'

def exponential(x, a, b, c):
    return a*np.exp(-x/b) + c

data_1 = []
data_2 = []
data_4 = []
for i in range(17):
    data = noiselib.loadmat(path.format(i))
    time = data['Qubit_Drive_to_Readout']
    data_1.append( data['Single_Shot_Occupation_SB1'] )
    data_2.append( data['Single_Shot_Occupation_SB2'] )
    data_4.append( data['Single_Shot_Occupation_SB3'] )
data_1 = np.vstack(data_1)
data_2 = np.vstack(data_2)
data_4 = np.vstack(data_4)
mean_1 = np.mean(data_1, axis=0)
mean_2 = np.mean(data_2, axis=0)
mean_4 = np.mean(data_4, axis=0)
popt_1, pcov_1 = curve_fit(exponential, time, mean_1, p0=[0.6, 40000, 0.2])
popt_2, pcov_2 = curve_fit(exponential, time, mean_2, p0=[0.6, 40000, 0.2])
popt_4, pcov_4 = curve_fit(exponential, time, mean_4, p0=[0.6, 40000, 0.2])

print(popt_1)
print(popt_2)
print(popt_4)

fig, ax = plt.subplots(1,1)
ax.plot(time, data_1.T, 'C0', alpha=0.3)
ax.plot(time, mean_1, label='Q1')
ax.plot(time, exponential(time, *popt_1), 'k:')
ax.plot(time, data_2.T, 'C1', alpha=0.3)
ax.plot(time, mean_2, label='Q2')
ax.plot(time, exponential(time, *popt_2), 'k:')
ax.plot(time, data_4.T, 'C2', alpha=0.3)
ax.plot(time, mean_4, label='Q4')
ax.plot(time, exponential(time, *popt_4), 'k:')
ax.legend()
plt.draw()
plt.pause(0.05)