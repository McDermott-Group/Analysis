import noiselib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np

def add_data_from_matlab(file_path):
    data = noiselib.loadmat(file_path)
    Amp = data['Amplitude']
    Charge = data['Q4_Charge_Bias']
    Freq = data['QB_Drive_Frequency']
    f_mean = np.mean(Freq)

    Amp = np.transpose(Amp)
    Charge = [float(np.format_float_positional(c, precision=2)) for c in Charge]
    Freq = [float(np.format_float_positional(f-f_mean, precision=4)) for f in Freq]

    return Amp, Charge, Freq
    # return np.log(Amp), Charge, Freq

file_path = "Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2" \
            "/06-27-21/Q4_spectroscopy_Charge/MATLABData/Q4_spectroscopy_Charge_002.mat"

Amp, Charge, Freq = add_data_from_matlab(file_path)
# print('Amp=', Amp)
# print('Charge=', Charge)
# print('Freq=', Freq)
plt.xticks(np.arange(0, 21, 1), Charge)
plt.yticks(np.arange(0, 17, 1)+0.35, Freq)
#
norm = LogNorm(vmin=np.min(Amp), vmax=np.max(Amp))
plt.imshow(Amp, cmap='binary', norm=norm)
# plt.imshow(Amp, cmap='binary')
plt.show()



