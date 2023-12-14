from junctionlib import *
from dataChest import *
from dataChest import dataChest as dc
import os

#DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data'
DATACHEST_ROOT = r'Z:\mcdermott-group\data'
base_path =  ['windowJJs', '20231120_ADR3']
device_name = '2000nm'

expt_paths =  [base_path + [device_name]]

paths = [os.path.join(*([DATACHEST_ROOT] + expt_path)) for expt_path in expt_paths]

for i,path in enumerate(paths):
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        # try:
        d = dataChest(os.path.join(path), dataSet)
        d.cd(expt_path)
        d.openDataset(dataSet,modify=False)
        junction = JJ(expt_path, [dataSet], [device_name])
        if "DIODE" in device_name.upper():
            junction.diodeCurrentFromOutputVoltage()
        #junction.plotLogIvsV(save=False)
        #junction.plotIvsV(save=False)
        junction.plotAppended(save=False) # added this
        junction.plotAppendedLogIvsV(save=False)
            # junction.autocenter()
        #except Exception as e:
            #print('Error Opening Dataset: {0}'.format(e))
            #continue
