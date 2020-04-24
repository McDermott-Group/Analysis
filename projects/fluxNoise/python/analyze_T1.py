import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# q, date, files = '\Q2', '\\03-24-20', range(251+1)
# q, date, files = '\Q2', '\\04-20-20', range(2,598+1)
# q, date, files = '\Q2', '\\04-20-20', range(641,973+1)
# q, date, files = '\Q2', '\\04-21-20', range(0,569+1)
# q, date, files = '\Q2', '\\04-23-20', range(4,244+1)
# q, date, files = '\Q1', '\\04-23-20', range(15,380+1)
# q, date, files = '\Q1', '\\04-24-20', range(0,73)
q, date, files = '\Q3', '\\04-24-20', range(0,655)
pathP1 = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar{}\General{}\P1\MATLABData\P1_{:03d}.mat'
pathT1 = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar{}\General{}\T1\MATLABData\T1_{:03d}.mat'
o = None
t = None
P1 = np.array([])

# for i in range(296+1):
for i in files:
    try:
        data = noiselib.loadmat( pathP1.format(q, date,i) )
    except:
        data = {'Single_Shot_Occupation': np.nan}
    P1 = np.append(P1, data['Single_Shot_Occupation'])
    
    data = noiselib.loadmat( pathT1.format(q, date,i) )
    o_new = np.array(data['Single_Shot_Occupation'])
    t = np.array(data['Qubit_Drive_to_Readout'])
    if o is None:
        o = np.empty((0,o_new.size))
    o = np.concatenate([o, [o_new]])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(files,P1)
ax.set_xlabel('File')
ax.set_ylabel('P1')
plt.draw()
plt.pause(0.05)

slice = np.mean(o[:,50:60], axis=1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(slice[P1 < 0.1],40)
ax.hist(slice[P1 > 0.1],40, alpha=0.6)
ax.set_xlabel('Single Shot Occupation near 70ns')
ax.set_ylabel('Counts')
plt.draw()
plt.pause(0.05)

thresh = 0.3 #0.24
above = o[slice > thresh].mean(axis=0)
below = o[slice < thresh].mean(axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Qubit Drive to Readout [ns]')
ax.set_ylabel('Single Shot Occupation')
ax.set_yscale('log')
# ax.plot(t, np.transpose(o), linewidth=0.2)
ax.plot(t, np.transpose(o[P1 > 0.1]), 'r', linewidth=0.2)
ax.plot(t, np.transpose(o[P1 < 0.1]), 'b', linewidth=0.2)
ax.plot(t, above, linewidth=4)
ax.plot(t, below, linewidth=4)
plt.draw()
plt.pause(0.05)

def exponential(t, b, a, tau):
    return b + a*np.exp(-t/tau)
    
def poisoned_exponential(t, b, a, tau_r, n_qp, tau_qp):
    return 0*b + a*np.exp(-t/tau_r)*np.exp(n_qp*(np.exp(-t/tau_qp)-1))

# popt_below_0, pcov_below_0 = curve_fit(exponential, t[:30], below[:30], p0=[0.1,1,30e3])
# popt_below_70, pcov_below_70 = curve_fit(exponential, t[-30:], below[-30:], p0=[0.1,1,30e3])
# popt_above_0, pcov_above_0 = curve_fit(exponential, t[:30], above[:30], p0=[0.1,1,30e3])
# popt_above_70, pcov_above_70 = curve_fit(exponential, t[-30:], above[-30:], p0=[0.1,1,30e3])

popt_below, pcov_below = curve_fit(poisoned_exponential, t, below, p0=[0.1,1.,30.e3,1.,100.e3])
popt_above, pcov_above = curve_fit(poisoned_exponential, t, above, p0=[0.1,1.,30.e3,1.,100.e3])

# ax.plot(t[:30], exponential(t[:30], *popt_below_0), 'k:', linewidth=2)
# ax.plot(t[-30:], exponential(t[-30:], *popt_below_70), 'k:', linewidth=2)
# ax.plot(t[:30], exponential(t[:30], *popt_above_0), 'k:', linewidth=2)
# ax.plot(t[-30:], exponential(t[-30:], *popt_above_70), 'k:', linewidth=2)
ax.plot(t, poisoned_exponential(t, *popt_below), 'k:', linewidth=2)
ax.plot(t, poisoned_exponential(t, *popt_above), 'k:', linewidth=2)
plt.draw()
plt.pause(0.05)

# print popt_below_0
# print popt_below_70
# print popt_above_0
# print popt_above_70

print popt_below
print popt_above