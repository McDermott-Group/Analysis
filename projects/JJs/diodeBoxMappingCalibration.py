from junctionlib import *
from dataChest import *
from dataChest import dataChest as dc
import os

# DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data'
DATACHEST_ROOT = 'Z:\\mcdermott-group\data'

base_path =  ['windowJJs', 'DiodeBoxMapping_IVData']

expt_paths =  [base_path]

paths = [os.path.join(*([DATACHEST_ROOT] + expt_path)) for expt_path in expt_paths]

reference_resistance = 172.9
for i,path in enumerate(paths):
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        try:
            print(path)
            print(dataSet)
            d = dc.dataChest(os.path.join(path), dataSet)
            d.cd(expt_path)
            d.openDataset(dataSet,modify=False)
            junction = JJ(expt_path, [dataSet], [''])
            junction._plot_I_vs_V_raw(save=False)

            if "Diode" not in dataSet:
                popt = junction._fit_reference_line(junction.V[0], junction.I[0])
                plt.figure()
                plt.title("Resistance = {:.2f}$\Omega$".format(1.0/popt[0]))
                plt.plot(junction.V[0],junction.I[0])
                plt.plot([-0.0003,0.0003],[popt[0]*-0.0003+popt[1],popt[0]*0.0003+popt[1]])
                plt.show()
                reference_resistance = popt[0]
                junction._diode_current_from_output_voltage()
            else:
                output_voltage = junction.I[0]*1e6
                output_current = junction.V[0]/reference_resistance
                plt.figure()
                plt.title("Diode Box Mapping")
                plt.xlabel("Output Current (A)")
                plt.ylabel("Output Voltage (V)")
                plt.plot(output_current,output_voltage)
                plt.show()
                np.savetxt('/Volumes/smb/mcdermott-group/data/windowJJs/Calibrations/20231201_cal.csv',np.vstack((output_voltage,output_current)).T,delimiter=',')
            # junction.autocenter()
        except Exception as e:
            print('Error Opening Dataset: {0}'.format(e))
            continue