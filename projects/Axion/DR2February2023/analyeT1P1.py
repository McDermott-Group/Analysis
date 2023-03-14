from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from axionlib import Parity,T1
import matplotlib
matplotlib.use('TkAgg')


base_path = ['Axion', '2023-02-26 - DR2']
user = 'DCH'
device_name = 'Axion4A Blackbody'
dates = ['03-03-23']
experiment_base_name ='T1_BB'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

for qb_id in [3,4,5,6]:
    t1s=[]
    p1s=[]
    temps=[]
    gammas=[]
    norm = [1e80,1e80,1e80,1e80,1e80,1e80,1e80]
    for i in range(len(paths)):
        path=paths[i]
        expt_path=expt_paths[i]
        for dataSet in os.listdir(path):
            try:
                d = dataChest(expt_path)
                d.openDataset(dataSet,modify=True)
                if d.getParameter('Qubit ID') != qb_id:
                    continue
                temps.append(d.getParameter('Blackbody Temperature')[0])
                print(dataSet)
                t1 = T1(expt_path, dataSet)
                fit_params = t1.fit(save=False)
                p1s.append(fit_params[2])
                # gammas.append((-1/fit_params[1])*fit_params[2]/(1-fit_params[2]))
                gammas.append((2e-6)*fit_params[2]/(1-fit_params[2]))
                if norm[qb_id] > (2e-6)*fit_params[2]/(1-fit_params[2]):
                    norm[qb_id] = (2e-6)*fit_params[2]/(1-fit_params[2])
                    print(norm)

                print(-1/fit_params[1])
            except:
                pass

    plt.semilogy([1000/temp for temp in temps],[norm[qb_id]/g for g in gammas],linestyle='',marker='*')
plt.xlabel('1/T (1/K)')
plt.ylabel('1 Hz/$\\Gamma_\\uparrow$ --Fixed T1')
plt.show()


