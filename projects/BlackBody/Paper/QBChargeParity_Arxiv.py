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
    Freq = [float(np.format_float_positional((f_mean-f)*1000, precision=4)) for f in Freq]

    return (Amp/np.min(Amp))**5, Charge, Freq
    # return np.log(Amp), Charge, Freq

file_path = "Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2" \
            "/06-27-21/Q4_spectroscopy_Charge/MATLABData/Q4_spectroscopy_Charge_002.mat"

Amp, Charge, Freq = add_data_from_matlab(file_path)
n_g = np.array([float(np.format_float_positional((c-0.4)*5/12, precision=2)) for c in Charge])

parity_x = np.arange(-5, 25, 0.1)
amp_p = 6
offsety = 8.35
offsetx = 2
period=12
parity_y_even = amp_p*np.cos((np.pi*(parity_x-offsetx)/period))+offsety
parity_y_odd = -amp_p*np.cos((np.pi*(parity_x-offsetx)/period))+offsety

plt.xticks(np.arange(2, 21, 6), n_g[2::6])
plt.yticks(np.arange(0, 17, 2)+0.49, Freq[::2])
plt.xlabel('$n_{g}$ (2e)')
plt.ylabel('$f-\overline{f_{01}}$ (MHz)')
#
# norm = LogNorm(vmin=np.min(Amp), vmax=np.max(Amp))
plt.plot(parity_x, parity_y_even, 'b')
plt.plot(parity_x, parity_y_odd, 'r')
plt.imshow(Amp, cmap='gray_r')
plt.show()



