import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
import TwoMeasDataFile
reload(TwoMeasDataFile)
from TwoMeasDataFile import *

"""This analyzes a single set of data to try to look at the white noise floor 
due to changes in idle time.  Every 10 datasets is a new idle time."""

filepath = ('Z:/mcdermott-group/data/fluxNoise2/DR1 - 2019-12-17/CorrFar/'
            'Q2/General/{}/QP_Tunneling_PSD/MATLABData/QP_Tunneling_PSD_{:03d}.mat')

# for j in range(6):
    # files = [filepath.format('06-11-20',i) for i in range(10*j,10*j+10)]
    # QPT = QPTunneling()
    # QPT.add_datasets(files, 'Single_Shot_Occupations_SB1')
    # QPT.plot_psd(figNum=111)
# plt.legend(['10','20','40','70','100','130'])

# for j in range(2):
    # files = [filepath.format('06-12-20',i) for i in range(2*j,2*j+2)]
    # QPT = QPTunneling()
    # QPT.add_datasets(files, 'Single_Shot_Occupations_SB1')
    # QPT.plot_psd(figNum=111)
    
    
filepath = ('Z:/mcdermott-group/data/fluxNoise2/DR1 - 2019-12-17/CorrFar/'
            'Q1Q2Corr/General/{}/Charge_resetting/MATLABData/Charge_resetting_{:03d}.mat')
files = [filepath.format('06-12-20',i) for i in [42,44,46,48,50,52]]
QPT = QPTunneling()
QPT.add_datasets(files, 'Single_Shot_Occupations_RO2_SB2')
QPT.plot_psd(figNum=111)