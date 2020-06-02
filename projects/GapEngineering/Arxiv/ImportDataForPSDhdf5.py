from dataChest import *
import matplotlib.pyplot as plt
# import os

file_path = 'GapEngineer\\Nb_GND_Dev01\\Leiden_2020Feb\\LIU\\Q1\\03-16-20\\QP_Tunneling_PSD\\HDF5Data'
file_name = 'cvd2133wum_QP_Tunneling_PSD.hdf5'

# ftag = file_name.split('_')[0]
print(file_path.split('\\'))
dc = dataChest(file_path.split('\\'))
dc.openDataset(file_name)

varsList = dc.getVariables()
data = dc.getData()
data = data.transpose()
pdata = data[4]
fig = plt.figure()
plt.plot(pdata, 'o-', label=r" Parity")
plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
plt.show()
