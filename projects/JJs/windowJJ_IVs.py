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
    for dataSet in os.listdir(path):
        # try:
        d = dataChest(os.path.join(path), dataSet)
        d.cd(expt_path)
        d.openDataset(dataSet,modify=False)
        junction = JJ(expt_path, [dataSet], [device_name])
        if "DIODE" in device_name.upper():
            junction._diode_current_from_output_voltage()
        junction.plot_I_vs_V(save=False, autocenter_mode='both',remove_jumps=True)
        junction.plot_I_vs_V(save=False, autocenter_mode='both', remove_jumps=False)
