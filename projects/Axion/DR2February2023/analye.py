from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from axionlib import Parity, T1
import matplotlib
matplotlib.use('TkAgg')

qb_id = 6

base_path = ['Axion', '2023-02-26 - DR2']
user = 'DCH'
device_name = 'Axion4A Blackbody'
dates = ['02-02-23']
experiment_base_name ='T1 BB'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

temperatures = []
P1s = []
T1s = []
for i in range(len(paths)):
    path = paths[i]
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        try:
            d = dataChest(os.path.join(path), dataSet)
            d.cd(expt_path)
            d.openDataset(dataSet,modify=True)
            if d.getParameter('Qubit ID') != qb_id:
                continue
        except:
            continue
        try:
            temperatures.append(d.getParameter('Blackbody Temperature'))
            T1_data = T1(path,dataSet)
        except:
            continue

plt.figure()
plt.semilogy([b*484*2*10 for b in biases],parity_rates,marker='.',linestyle='None')
plt.xlabel('Radiator Frequency')
plt.ylabel('Parity Rate (GHz)')
plt.show()


