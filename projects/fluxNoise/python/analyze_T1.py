import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# q, date, files, thresh = '\Q2', '\\03-24-20', range(251+1)
# q, date, files, thresh = '\Q2', '\\04-20-20', range(2,598+1)
# q, date, files, thresh = '\Q2', '\\04-20-20', range(641,973+1)
# q, date, files, thresh = '\Q2', '\\04-21-20', range(0,569+1), 0.31 ###
# q, date, files, thresh = '\Q2', '\\04-23-20', range(4,244+1), 0.375
# q, date, files, thresh = '\Q1', '\\04-23-20', range(15,380+1), 0.1
# q, date, files = '\Q1', '\\04-24-20', range(0,73)
# q, date, files, thresh = '\Q3', '\\04-24-20', range(400,999), 0.026
# q, date, files, thresh = '\Q4', '\\04-24-20', range(0,357), 0.07
q, date, files, thresh = '\Q4', '\\04-25-20', range(0,641), 0.19
# q, date, files, thresh = '\Q4', '\\04-26-20', range(0,835), 0.26
pathP1 = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar{}\General{}\P1\MATLABData\P1_{:03d}.mat'
pathT1 = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar{}\General{}\T1\MATLABData\T1_{:03d}.mat'
o = None
t = None
P1 = np.array([])

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

slice_index = 35
slice = np.mean(o[:,slice_index-5:slice_index+5], axis=1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(files,P1, linewidth=0.5)
ax.scatter(files,P1,c=slice>thresh, cmap='bwr')
ax.set_xlabel('File')
ax.set_ylabel('P1')
plt.draw()
plt.pause(0.05)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(slice,40)
# ax.hist(slice[P1 < thresh],40)
# ax.hist(slice[P1 > thresh],40, alpha=0.6)
ax.set_xlabel('Single Shot Occupation near {:.2f}us'.format(t[slice_index]/1000))
ax.set_ylabel('Counts')
plt.draw()
plt.pause(0.05)

above = o[slice > thresh].mean(axis=0)
below = o[slice < thresh].mean(axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Qubit Drive to Readout [ns]')
ax.set_ylabel('Single Shot Occupation')
ax.set_yscale('log')
# ax.plot(t, np.transpose(o), linewidth=0.2)
ax.plot(t, np.transpose(o[slice > thresh]), 'r', linewidth=0.2)
ax.plot(t, np.transpose(o[slice < thresh]), 'b', linewidth=0.2)
ax.plot(t, above, linewidth=4)
ax.plot(t, below, linewidth=4)
plt.draw()
plt.pause(0.05)

def exponential(t, b, a, tau):
    return b + a*np.exp(-t/tau)
    
def poisoned_exponential(t, b, a, tau_r, n_qp, tau_qp):
    return b + a*np.exp(-t/tau_r)*np.exp(n_qp*(np.exp(-t/tau_qp)-1))
    
def poisoned_exponential_simul(t, a, tau_r, n_qp0, n_qp1, tau_qp0, tau_qp1):
    return np.ravel([poisoned_exponential(t, 0, a, tau_r, n_qp0, tau_qp0),
                     poisoned_exponential(t, 0, a, tau_r, n_qp1, tau_qp1)])

# popt_below, pcov_below = curve_fit(poisoned_exponential, t, below, p0=[0.,1.,30.e3,1.,100.e3])
# popt_above, pcov_above = curve_fit(poisoned_exponential, t, above, p0=[0.,1.,30.e3,1.,100.e3])
popt_simul, pcov_simul = curve_fit(poisoned_exponential_simul, 
                                   t, np.ravel([below,above]), 
                                   p0=[1.,30.e3,1.,1.,100.e3,100.e3])

fits = poisoned_exponential_simul(t, *popt_simul).reshape((2,t.size))
ax.plot(t, np.transpose(fits), 'k:', linewidth=2)
plt.draw()
plt.pause(0.05)

print popt_simul