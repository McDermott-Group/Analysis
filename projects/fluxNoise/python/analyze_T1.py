import os
import numpy as np
import noiselib
reload(noiselib)
from noiselib import matpaths
from dataChest import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
# q, date, filesP1_I, filesP1_X, filesT1, thresh = 'Q3', '04-29-20', range(112,493), range(112,493), range(0,381), 0.33
# q, date, filesP1_I, filesP1_X, filesT1, thresh = 'Q3', '04-29-20', range(0,63), range(0,63), range(0,63), 0.25
q, date, filesP1_I, filesP1_X, filesT1, thresh = 'Q4', '04-30-20', range(0,231), range(0,231), range(0,231), 0.44

fP1_I = matpaths(q,date,'P1_I',filesP1_I)
fP1_X = matpaths(q,date,'P1_X',filesP1_X)
fT1 = matpaths(q,date,'T1',filesT1)
o = None
t = None
P1_I = np.array([])
P1_X = np.array([])

for f in fP1_I:
    try:
        data = noiselib.loadmat( f )
    except:
        data = {'Single_Shot_Occupation': np.nan}
    P1_I = np.append(P1_I, data['Single_Shot_Occupation'])
for f in fP1_X:
    try:
        data = noiselib.loadmat( f )
    except:
        data = {'Single_Shot_Occupation': np.nan}
    P1_X = np.append(P1_X, data['Single_Shot_Occupation'])
for f in fT1:
    data = noiselib.loadmat( f )
    o_new = np.array(data['Single_Shot_Occupation'])
    t = np.array(data['Qubit_Drive_to_Readout'])
    if o is None:
        o = np.empty((0,o_new.size))
    o = np.concatenate([o, [o_new]])

slice_index = 45
slice = np.mean(o[:,slice_index-5:slice_index+5], axis=1)
# filter = slice > thresh
filter = P1_I > 0.064

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(filesP1_I, P1_I, linewidth=0.5)
ax.scatter(filesP1_I, P1_I, c=filter, cmap='bwr', s=4)
ax.plot(filesP1_X, P1_X, linewidth=0.5)
ax.scatter(filesP1_X, P1_X, c=filter[:len(P1_X)], cmap='bwr', s=4)
ax.set_xlabel('File')
ax.set_ylabel('P1')
plt.draw()
plt.pause(0.05)

fig = plt.figure()
ax = fig.add_subplot(111)
h, bins = np.histogram(slice, bins=40)
ax.hist(slice[filter], bins, color='red')
ax.hist(slice[~filter], bins, color='blue', alpha=0.6)
ax.set_xlabel('Single Shot Occupation near {:.2f}us'.format(t[slice_index]/1000))
ax.set_ylabel('Counts')
plt.draw()
plt.pause(0.05)

yfilter = o[filter].mean(axis=0)
nfilter = o[~filter].mean(axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Qubit Drive to Readout [ns]')
ax.set_ylabel('Single Shot Occupation')
ax.set_yscale('log')
# ax.plot(t, np.transpose(o), linewidth=0.2)
ax.plot(t, np.transpose(o[filter]), 'r', linewidth=0.2)
ax.plot(t, np.transpose(o[~filter]), 'b', linewidth=0.2)
ax.plot(t, yfilter, linewidth=4)
ax.plot(t, nfilter, linewidth=4)
plt.draw()
plt.pause(0.05)

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

print '{:>10}{:>10}{:>10}'.format('','red','blue')
print '{:>10}{:10.3f}{:10.3f}'.format('A',*popt[0:2])
print '{:>10}{:10.3f}{:10}'.format('tau_R',popt[2]*1e-3,'')
print '{:>10}{:10.3f}{:10.3f}'.format('n_QP',*popt[3:5])
print '{:>10}{:10.3f}{:10.3f}'.format('tau_QP',popt[5]*1e-3,popt[6]*1e-3)