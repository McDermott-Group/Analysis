from IQfit_QM import *
from dataChest import *
import matplotlib
matplotlib.use('TkAgg')

#path = r"Z:\mcdermott-group\data\Micromachining\2023-01-23 - DR2\MAS\MM_SC_1\02-06-23\Resonator_Spectroscopy"
path = ['Micromachining', '2023-01-23 - DR2', 'MAS','MM_SC_1','02-07-23', 'Resonator_Spectroscopy']

#IQ fit of file and frequency guess
Iq1 = IQ_fit(path, 6.427)

print(Iq1.res_freq)
Iq1.IQ_fit_single_power(file_index=2, power=10, Qc_guess=None, Qi_guess=None, Singleshotlike=True, ind=1, auto_res_freq=True)
# Iq1.IQ_fit_multi_power(powers=np.linspace(20,-12.5,14), Qc_guess = 5e6, att=40, Qi_guess=None, save=True, ind=1, auto_res_freq=False)
# print(os.path.basename(p))


# p = r"Axion\2022-01-04 - ADR3\Dave\Axion1A Qubit\01-03-22\CavitySpectroscopy_Power\HDF5Data"
#12,20 - 6.6 - 4.7k - A+
#18 - 6.7 - 5.25k - A+
#22 - 6.5 - 2.4k - A+

#26 - 6.25 - 5.9k - A+
#32,38 - 6.15 - 4k - B
#40 - 6.05 - 3.4k - C


'''
#Z:\mcdermott-group\data\Micromachining\2022-09-22 - DR2\MAS\MM_1\11-10-22\T1_Poisoning_vs_bias_vary_idle
path = ['Micromachining', '2023-01-23 - DR2', 'MAS','MM_SC_1','02-06-23', 'Resonator_Spectroscopy']
dc = dataChest(path)
contents = dc.ls()
files = contents[0]
f = files[0]
print(f)
dc.openDataset(f)
varsList = dc.getVariables()
print(varsList)
'''