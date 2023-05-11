from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from axionlib import Parity
import matplotlib
from tqdm import tqdm
matplotlib.use('TkAgg')

qb_id = 4
rad_id = 1

base_path = ['SW_Circmon', '2023-02-15 - HQAN']
user = 'SW'
device_name = 'SW_Circmon'
dates = ['02-17-23']
experiment_base_name ='Qubit_Spectroscopy_Poormans_PSD'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

f01s = []
for i in tqdm(range(len(paths))):
    path = paths[i]
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        try:
            d = dataChest(os.path.join(path), dataSet)
            d.cd(expt_path)
            d.openDataset(dataSet, modify=False)
        except Exception as e:
            print(e)
        try:
            data=d.getData()
            f01s.append(d.getParameter("Q1 Peak frequency")[0])
            print(d.getParameter("Q1 Peak frequency")[0])
            # idx=np.argmax(data[:,3])
            # f01s.append(data[idx][0])
            # print(data[idx][0])
        except Exception as e:
            print(e)

plt.figure()
plt.hist(f01s,bins='auto')
plt.show()
plt.figure()
plt.plot(f01s, ".")
plt.show()