import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q2\General\03-24-20\T1\MATLABData\T1_{:03d}.mat'
o = None
t = None

# for i in range(296+1):
for i in range(251+1):
    data = noiselib.loadmat( path.format(i) )
    o_new = np.array(data['Single_Shot_Occupation'])
    t = np.array(data['Qubit_Drive_to_Readout'])
    if o is None:
        o = np.empty((0,o_new.size))
    o = np.concatenate([o, [o_new]])

slice = np.mean(o[:,65:75], axis=1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(slice,40)
ax.set_xlabel('Single Shot Occupation near 70ns')
ax.set_ylabel('Counts')
plt.draw()
plt.pause(0.05)

thresh = 0.22
above = o[slice > thresh].mean(axis=0)
below = o[slice < thresh].mean(axis=0)

def exponential(t, b, a, tau):
    return b + a*np.exp(-t/tau)

popt_below_0, pcov_below_0 = curve_fit(exponential, t[:30], below[:30], p0=[0.1,1,30e3])
popt_below_70, pcov_below_70 = curve_fit(exponential, t[-30:], below[-30:], p0=[0.1,1,30e3])
popt_above_0, pcov_above_0 = curve_fit(exponential, t[:30], above[:30], p0=[0.1,1,30e3])
popt_above_70, pcov_above_70 = curve_fit(exponential, t[-30:], above[-30:], p0=[0.1,1,30e3])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, np.transpose(o), linewidth=0.2)
ax.plot(t, above, linewidth=4)
ax.plot(t, below, linewidth=4)
ax.plot(t[:30], exponential(t[:30], *popt_below_0), 'k:', linewidth=2)
ax.plot(t[-30:], exponential(t[-30:], *popt_below_70), 'k:', linewidth=2)
ax.plot(t[:30], exponential(t[:30], *popt_above_0), 'k:', linewidth=2)
ax.plot(t[-30:], exponential(t[-30:], *popt_above_70), 'k:', linewidth=2)
ax.set_xlabel('Qubit Drive to Readout [ns]')
ax.set_ylabel('Single Shot Occupation')
ax.set_yscale('log')
plt.draw()
plt.pause(0.05)

print popt_below_0
print popt_below_70
print popt_above_0
print popt_above_70