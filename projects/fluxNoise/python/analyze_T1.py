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

# q, date, files, thresh = 'Q2', '03-24-20', range(251+1)
# q, date, files, thresh = 'Q2', '04-20-20', range(2,598+1)
# q, date, files, thresh = 'Q2', '04-20-20', range(641,973+1)
# q, date, files, thresh = 'Q2', '04-21-20', range(0,569+1), 0.31 ###
# q, date, files, thresh = 'Q2', '04-23-20', range(4,244+1), 0.375
# q, date, files, thresh = 'Q1', '04-23-20', range(15,380+1), 0.1
# q, date, files = 'Q1', '04-24-20', range(0,73)
# q, date, files, thresh = 'Q3', '04-24-20', range(400,999), 0.026
# q, date, files, thresh = 'Q4', '04-24-20', range(0,357), 0.07
# q, date, files, thresh = 'Q4', '04-25-20', range(0,641), 0.19 ###
# q, date, files, thresh = 'Q4', '04-26-20', range(0,835), 0.26
# q, date, filesP1, filesT1, thresh = 'Q3', '04-29-20', range(112,493), range(0,381), 0.33 # really bad
# q, date, filesP1, filesT1, thresh = 'Q3', '04-29-20', range(0,63), range(0,63), 0.25 # shows nothing

# ds = q4_0430_T1
# ds = q3_0501_T1
# ds = q1_0503_T1
ds = q2_0505_T1

fP1_I = matpaths(fileName='P1_I', fileNums='files_P1', **ds)
fP1_X = matpaths(fileName='P1_X', fileNums='files_P1', **ds)
fT1 = matpaths(fileName='T1', fileNums='files_T1', **ds)
o = None
t = None
P1 = { 'I':np.array([]),
       'X': np.array([]) }

if 'recalibrate' in ds and ds['recalibrate']:
    print 'recalibrating...'
    data = noiselib.loadmat( fP1_I[0] )
    iq0 = np.array((data['Is'],data['Qs']))
    data = noiselib.loadmat( fP1_X[0] )
    iq1 = np.array((data['Is'],data['Qs']))
    c = cal.CalibratedStates( (0,iq0), (1,iq1), plot=False)

for files, gate in [ (fP1_I,'I'), (fP1_X,'X') ]:
    for f in files:
        try:
            data = noiselib.loadmat( f )
        except:
            print 'corrupted file:', f
            data = { 'Single_Shot_Occupations': [np.nan],
                     'Single_Shot_Occupation': [np.nan]}
        if 'recalibrate' in ds and ds['recalibrate']:
            states, SSO, _ = c.get_single_shot_occupation( np.array((data['Is'],
                                                                     data['Qs'])) )
            sso = SSO[1]
        else:
            sso = data['Single_Shot_Occupation']
        P1[gate] = np.append(P1[gate], sso)
for f in fT1:
    data = noiselib.loadmat( f )
    o_new = np.array(data['Single_Shot_Occupation'])
    t = np.array(data['Qubit_Drive_to_Readout'])
    if o is None:
        o = np.empty((0,o_new.size))
    o = np.concatenate([o, [o_new]])

# take a slice of all the T1 curves, filter based on this or on P1
slice_index = 45
slice = np.mean(o[:,slice_index-5:slice_index+5], axis=1)
if 'thresh_P1' in ds:
    filter = P1['I'] > ds['thresh_P1']
elif 'thresh_T1' in ds:
    filter = slice > ds['thresh_T1']

# plot P1
fig, ax = plt.subplots()
x = range(len(fP1_I))
ax.plot(x, P1['I'], linewidth=0.5)
ax.scatter(x, P1['I'], c=filter, cmap='bwr', s=4)
ax.plot(x, P1['X'], linewidth=0.5)
ax.scatter(x, P1['X'], c=filter[:len(P1['X'])], cmap='bwr', s=4)
ax.set_xlabel('File')
ax.set_ylabel('P1')
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
plt.draw()
plt.pause(0.05)

# hist of T1 cut
fig, ax = plt.subplots()
h, bins = np.histogram(slice, bins=40)
ax.hist(slice[filter], bins, color='red')
ax.hist(slice[~filter], bins, color='blue', alpha=0.6)
ax.set_xlabel('Single Shot Occupation near {:.2f}us'.format(t[slice_index]/1000))
ax.set_ylabel('Counts')
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
plt.draw()
plt.pause(0.05)

yfilter = o[filter].mean(axis=0)
nfilter = o[~filter].mean(axis=0)

# T1 curves
fig, ax = plt.subplots()
ax.set_xlabel('Qubit Drive to Readout [ns]')
ax.set_ylabel('Single Shot Occupation')
ax.set_yscale('log')
fig.suptitle(ds['Q']+'\n'+str(ds['date']))
# ax.plot(t, np.transpose(o), linewidth=0.2)
ax.plot(t, np.transpose(o[filter]), 'r', linewidth=0.2)
ax.plot(t, np.transpose(o[~filter]), 'b', linewidth=0.2)
ax.plot(t, yfilter, linewidth=4)
ax.plot(t, nfilter, linewidth=4)
plt.draw()
plt.pause(0.05)

# fit T1 and plot fits

def exponential(t, b, a, tau):
    return b + a*np.exp(-t/tau)
    
def poisoned_exponential(t, b, a, tau_r, n_qp, tau_qp):
    return b + a*np.exp(-t/tau_r)*np.exp(n_qp*(np.exp(-t/tau_qp)-1))
    
def poisoned_exponential_simul(t, a0, a1, tau_r, n_qp0, n_qp1, tau_qp0, tau_qp1):
    return np.ravel([poisoned_exponential(t, 0, a0, tau_r, n_qp0, tau_qp0),
                     poisoned_exponential(t, 0, a1, tau_r, n_qp1, tau_qp1)])

popt, pcov = curve_fit(poisoned_exponential_simul, 
                       t, np.ravel([yfilter,nfilter]), 
                       p0=[1.,1.,30.e3,1.,1.,100.e3,100.e3])

fits = poisoned_exponential_simul(t, *popt).reshape((2,t.size))
ax.plot(t, np.transpose(fits), 'k:', linewidth=2)
plt.draw()
plt.pause(0.05)

text  = '{:>7}{:>10}{:>10}\n'.format('','red','blue')
# text += '{:>7}{:10.3f}{:10.3f}\n'.format('A',*popt[0:2])
text += '{:>7}{:>10.3f}{:>10}\n'.format('tau_R',popt[2]*1e-3,'')
text += '{:>7}{:>10.3f}{:>10.3f}\n'.format('n_QP',*popt[3:5])
text += '{:>7}{:>10.3f}{:>10.3f}\n'.format('tau_QP',popt[5]*1e-3,popt[6]*1e-3)
ax.text(0, o[o>0].min(), text, family='monospace')
plt.draw()
plt.pause(0.05)