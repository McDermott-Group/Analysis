import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
import matplotlib.pyplot as plt
import noiselib
from numba import jit
import time
import impact_lib
import importlib
importlib.reload(impact_lib)
from impact_lib import *
import pickle
import gc

path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'

if __name__ == '__main__':
    app = Controller(sys.argv, fQ=1., calc_all=False, plot=True,
                     pdfs_file=['{}/ChargePDFs_{}.npy'.format(path,300),
                                '{}/ChargePDFs_{}.npy'.format(path,100)],
                     event_files=[path+'/Gamma.txt'])#, 
                                  # path+'/Gamma_10deg.txt'],
                                  # path+'/Gamma_10deg_pt2.txt'])
    
    if hasattr(app, 'view'):
        sys.exit(app.exec_())