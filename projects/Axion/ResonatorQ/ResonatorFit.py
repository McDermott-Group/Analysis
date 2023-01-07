from IQfit import *


p = r"Axion\2022-09-08 - ADR3\Dave\Axion3A_20220802D_Ox3\09-08-22\CavitySpectroscopy_Power\HDF5Data"

Iq1 = IQ_fit(p, 6.809)
print(Iq1.res_freq)
Iq1.IQ_fit_single_power(file_index=8, power=10, Qc_guess=None, Qi_guess=None, CavitySpectroscopylike=True, ind=1, auto_res_freq=True)
# Iq1.IQ_fit_multi_power(powers=np.linspace(20,-12.5,14), Qc_guess = 5e6, att=40, Qi_guess=None, save=True, ind=1, auto_res_freq=False)
# print(os.path.basename(p))


# p = r"Axion\2022-01-04 - ADR3\Dave\Axion1A Qubit\01-03-22\CavitySpectroscopy_Power\HDF5Data"
#12,20 - 6.6 - 4.7k - A+
#18 - 6.7 - 5.25k - A+
#22 - 6.5 - 2.4k - A+

#26 - 6.25 - 5.9k - A+
#32,38 - 6.15 - 4k - B
#40 - 6.05 - 3.4k - C