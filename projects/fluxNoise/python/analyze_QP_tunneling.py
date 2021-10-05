import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean
import QPTunneling
reload(QPTunneling)
from QPTunneling import *
import ChargeOffset
reload(ChargeOffset)
from ChargeOffset import *
import TwoMeasDataFile
reload(TwoMeasDataFile)
from TwoMeasDataFile import *

""""""

filepath = ('Z:/mcdermott-group/data/fluxNoise/DR1 - 2019-08-12/{}/'
            '{}/General/{}/QP_Tunneling_PSD/MATLABData/QP_Tunneling_PSD_{:03d}.mat')

data = [
('Circ2', 'Q1', '08-26-19', np.arange(0,1113+1)),
('Circ2', 'Q3', '08-26-19', np.arange(0,2536+1)),
('Circ1', 'Q2', '08-30-19', np.arange(0,1242+1)),
('Circ2', 'Q4', '09-07-19', np.arange(0,264+1)),
('Circ1', 'Q3', '09-16-19', np.arange(0,263+1))
]

fig = plt.figure()
for d in data:
    QPT = QPTunneling()
    files = [filepath.format(d[0],d[1],d[2],i) for i in d[3]]
    QPT.add_datasets(files)
    QPT.plot_psd(figNum=fig.number, label=d[0]+'-'+d[1])