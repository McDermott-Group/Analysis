import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
from dataChest import dataChest
from random import randrange


date = '08-13-20'
QP_Neg_4dBm_path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/Outside/LIU/Q4_withQ5Poison/{}/QP_Tunneling_PSD_Poison_Negatvie_4dBm/MATLABData/')
QP_Neg_4dBm_files = np.arange(0,19, 1) # 299 is good
QP_Neg_4dBm_filenames = [QP_Neg_4dBm_path.format(date) + 'QP_Tunneling_PSD_Poison_Negatvie_4dBm_{:03d}.mat'.format(i) for i in QP_Neg_4dBm_files]
QP_Neg_4dBm_psd = QP_PSD()
QP_Neg_4dBm_psd.add_data_from_matlab(QP_Neg_4dBm_filenames)
QP_Neg_4dBm_psd.plot_PSD(fit=False)

