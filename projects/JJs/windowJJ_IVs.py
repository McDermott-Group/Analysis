from junctionlib import *
from dataChest import *
import os
import platform

DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data/' if 'macOS' in platform.platform() else 'Z:\\mcdermott-group\\data\\'


base_path =  ['windowJJs', '20231208_ADR3']
device_name = '2000nm (Bottom Left)'
# '325nm (Top Right)'  || '350nm (Top Left)' || '350nm (Top Right)' || '2000nm (Bottom Left)'
# '2000nm (Top Left)'
expt_paths =  [base_path + [device_name]]

paths = [os.path.join(*([DATACHEST_ROOT] + expt_path)) for expt_path in expt_paths]

for i,path in enumerate(paths):
    expt_path = expt_paths[i]
    for data_set in os.listdir(path):
        try:
            d = dataChest(os.path.join(path), data_set)
            d.cd(expt_path)
            d.openDataset(data_set, modify=False)
            junction = JJ(expt_path, [data_set], [device_name])
            if "DIODE" in device_name.upper():
                junction.diodeCurrentFromOutputVoltage()
            junction.plotLogIvsV(save=False,autocenter_mode='gap')
            # junction.plotIvsV(save=False,autocenter_mode='gap')
        except Exception as e:
            print('Error Opening Dataset: {0}'.format(e))
            continue
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
