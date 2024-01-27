from datetime import datetime
from scipy import signal
from dataChest import dataChest as dc
import os
import matplotlib.pyplot as plt
from axionlib import *
from scipy.optimize import curve_fit
import numpy as np

import matplotlib
matplotlib.use('TkAgg')

qb_id = 3
rad_id = None

data_root = r'/Volumes/smb/mcdermott-group/data'#r'C:\Local Data'

base_path =  ['Axion', '2024-01-07 - DR5']#['Axion', '2023-04-05 - DR2']
user = 'DCH'
device_name = 'Axion4A DieC'#'Axion4A Blackbody'
dates = ['01-10-24']#['04-16-23']
experiment_base_name ='PSD_Q3_Fast'#'PSD_Q2_99mK_BB'#'PSD_Q{:d}_J{:d}B'.format(qb_id,rad_id)

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]

paths = [os.path.join(*([data_root] + expt_path)) for expt_path in expt_paths]

time_per_iteration = 100000#10000000

psds = []
fit_psds = []

parity_rates=[]
fidelities=[]
biases=[]
for i in range(len(paths)):
    path = paths[i]
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        try:
            d = dc.dataChest(os.path.join(path), dataSet)
            d.cd(expt_path)
            d.openDataset(dataSet,modify=True)
        except Exception as e:
            print(e)
            print('Error Opening Dataset')
            continue
        try:
            params = d.getParameterList()
            if 'Radiator ID' not in params:
                d.addParameter('Radiator ID',0)
                d.addParameter('JB0 Voltage Bias',0,'mV')
            if rad_id is not None:
                if 'JB{:1d} Voltage Bias'.format(rad_id)[0] not in params:
                    biases.append(d.getParameter('JB{:1d} Voltage Bias'.format(rad_id))[0])
                else:
                    biases.append(0)
            else:
                biases.append(0)
            if 'Qubit ID' not in params:
                d.addParameter('Qubit ID',qb_id)
            if 'Time Per Iteration' not in params:
                d.addParameter('Time Per Iteration', time_per_iteration,'ns',overwrite=False)
            try:
                # raise Exception('re-analyze please')
                parity_rates.append(d.getParameter('Fit Parity Rate')[0])
                fidelities.append(d.getParameter('Fit Fidelity'))
            except:
                p=Parity(expt_path, dataSet, rad_id)
                f,psd = p.psd()
                psds.append(psd)
                f, fit_psd = p.fit_psd()
                parity_rate,fidelity = p.fit()
                parity_rates.append(parity_rate)
                fidelities.append(fidelity)
                p.plot(fit=True,save=True)
        except Exception as e:
            print(e)
            continue

#
# p._psd = np.average(np.array(psds),axis=0)
# p._fit_psd = None
# p.plot()
#
# plt.show()
#
# plt.figure()
# plt.semilogy([(b)*484*2*10 for b in biases],parity_rates,marker='.',linestyle='None',color='b')
# plt.gca().xaxis.set_tick_params(direction = 'in')
# plt.gca().yaxis.set_tick_params(direction = 'in')
# plt.xlabel('Radiator Frequency (GHz)')
# plt.ylabel('Parity Switching Rate (GHz)')
# plt.show()
#
# filtered_rates=[]
# filtered_biases=[]
# for i,rate in enumerate(parity_rates):
#     if rate<1500:
#         filtered_rates.append(rate)
#         filtered_biases.append(biases[i])
#
# plt.figure()
# plt.semilogy([(b)*484*2*10 for b in filtered_biases],filtered_rates,marker='.',linestyle='None',color='b')
# plt.gca().xaxis.set_tick_params(direction = 'in')
# plt.gca().yaxis.set_tick_params(direction = 'in')
# plt.xlabel('Radiator Frequency (GHz)')
# plt.ylabel('Parity Switching Rate (GHz)')
# plt.show()
#
#
#

