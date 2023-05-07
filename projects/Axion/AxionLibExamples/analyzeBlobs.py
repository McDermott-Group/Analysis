from axionlib import *
from dataChest import *
import os

base_path = ['Axion', '2023-02-26 - DR2']
user = 'DCH'
device_name = 'Axion4A Blackbody'
dates = ['03-02-23']
experiment_base_name ='P1'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

qb_data = {i:[] for i in range(1,7)}
for path,expt_path in zip(paths,expt_paths):
    for dataSet in os.listdir(path):
        d = dataChest(expt_path)
        d.openDataset(dataSet, modify=True)
        qb_data[d.getParameter('Qubit ID')].append(dataSet)

for qb in range(1,7):
    for expt_path in expt_paths:
        try:
            IQBlobs(expt_path,qb_data[qb]).plot(save=True)
        except:
            continue