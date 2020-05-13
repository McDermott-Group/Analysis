import os
import numpy as np
import noiselib
reload(noiselib)
from noiselib import matpaths
from dataChest import *
import general.calibration as cal
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datasets
reload(datasets)
from datasets import *

# q, date, filesP1, filesT1, thresh = 'Q3', '04-29-20', range(112,493), range(0,381), 0.33 # really bad
# q, date, filesP1, filesT1, thresh = 'Q3', '04-29-20', range(0,63), range(0,63), 0.25 # shows nothing
# q, date, filesP1, filesT1, thresh = 'Q4', '04-30-20', range(0,231), range(0,231), 0.44

#################################

# ds = q4_0430_T1
# ds = q3_0501_T1
# ds = q1_0503_T1
ds = q2_0505_T1

fP1_I = matpaths(fileName='P1_I', fileNums='files_P1', **ds)
fP1_X = matpaths(fileName='P1_X', fileNums='files_P1', **ds)
dwellsDirty = np.array([])
dwellsClean = np.array([])

dwellsDirty1 = np.array([])
dwellsDirty2 = np.array([])
dwellsDirty3 = np.array([])

if 'recalibrate' in ds and ds['recalibrate']:
    data = noiselib.loadmat( fP1_I[0] )
    iq0 = np.array((data['Is'],data['Qs']))
    data = noiselib.loadmat( fP1_X[0] )
    iq1 = np.array((data['Is'],data['Qs']))
    c = cal.CalibratedStates( (0,iq0), (1,iq1), plot=False)

for i,f in enumerate(fP1_I):
    try:
        data = noiselib.loadmat( f )
    except:
        print 'corrupted file:', f
        data = {'Single_Shot_Occupations': [np.nan],
                'Single_Shot_Occupation': [np.nan]}
    if 'recalibrate' in ds and ds['recalibrate']:
        states, SSO, _ = c.get_single_shot_occupation( np.array((data['Is'],
                                                                 data['Qs'])) )
    else:
        states = data['Single_Shot_Occupations']
        SSO = {1: data['Single_Shot_Occupation']}
    d0,d1 = noiselib.bin_dwell_times(states)
    if SSO[1] > ds['thresh_P1']:
        dwellsDirty = np.append(dwellsDirty, d0)
        if i in range(30,40):
            dwellsDirty1 = np.append(dwellsDirty, d0)
        elif i in range(80,88):
            dwellsDirty2 = np.append(dwellsDirty, d0)
        elif i in range(311,319):
            dwellsDirty3 = np.append(dwellsDirty, d0)
    else:
        dwellsClean = np.append(dwellsClean, d0)

fig, ax = plt.subplots()
ax.hist(dwellsDirty1, bins, color='red')
ax.hist(dwellsDirty2, bins, color='blue', alpha=0.6)
ax.hist(dwellsDirty3, bins, color='green', alpha=0.6)
ax.set_xlabel('Dwell Time in 0 State Before Excitation [Reps=500us]')
ax.set_ylabel('Counts')
ax.set_yscale('log')
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
plt.draw()
plt.pause(0.05)

fig, ax = plt.subplots()
hD, bins, _ = ax.hist(dwellsDirty,30, color='red')
hC, bins, _ = ax.hist(dwellsClean, bins, color='blue', alpha=0.6)
ax.set_xlabel('Dwell Time in 0 State Before Excitation [Reps=500us]')
ax.set_ylabel('Counts')
ax.set_yscale('log')
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
plt.draw()
plt.pause(0.05)

def exponential(t, a, tau):
    return a*np.exp(-t/tau)
    
def poisoned_exponential(t, b, a, tau_r, n_qp, tau_qp):
    return b + a*np.exp(-t/tau_r)*np.exp(n_qp*(np.exp(-t/tau_qp)-1))
    
def poisoned_exponential_simul(t, a0, a1, tau_r, n_qp0, n_qp1, tau_qp0, tau_qp1):
    return np.ravel([poisoned_exponential(t, 0, a0, tau_r, n_qp0, tau_qp0),
                     poisoned_exponential(t, 0, a1, tau_r, n_qp1, tau_qp1)])

t = (bins[1:]+bins[:-1])/2

poptD, pcovD = curve_fit(exponential, t[t.size/2:], hD[t.size/2:], p0=[hD.max()/2, .3])
poptC, pcovC = curve_fit(exponential, t[t.size/2:], hC[t.size/2:], p0=[hC.max()/2, .3])

ax.plot(t, exponential(t, *poptD), 'k:')
ax.plot(t, exponential(t, *poptC), 'k:')
ax.set_ylim([0.5, None])
plt.draw()
plt.pause(0.05)

popt, pcov = curve_fit(poisoned_exponential_simul, 
                       t, np.ravel([hD,hC]), 
                       p0=[poptD[0],poptC[0],poptC[1],1.,0.,.5,.5], 
                       bounds=(0,np.inf))
                       # bounds=([0,0,0,0,0,0,0],[np.inf,np.inf,np.inf,3,3,np.inf,np.inf]))

fits = poisoned_exponential_simul(t, *popt).reshape((2,t.size))
ax.plot(t, np.transpose(fits), 'k', linewidth=2)
plt.draw()
plt.pause(0.05)

text  = '{:>7}{:>12}{:>12}\n'.format('','Dirty','Clean')
# text += '{:>7}{:>12.3f}{:>12.3f}\n'.format('A',*popt[0:2])
text += '{:>7}{:>12.3f}{:>12}\n'.format('tau_R',popt[2]*500,'')
text += '{:>7}{:>12.3f}{:>12.3f}\n'.format('n_QP',*popt[3:5])
text += '{:>7}{:>12.3f}{:>12.3f}\n'.format('tau_QP',popt[5]*500,popt[6]*500)
text += '{:>7}{:>12.3f}{:>12.3f}\n'.format('tau',poptD[1]*500,poptC[1]*500)
ax.text(t[-1],fits.max()/10, text, family='monospace', horizontalalignment='right')
plt.draw()
plt.pause(0.05)