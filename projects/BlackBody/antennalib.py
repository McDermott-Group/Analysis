"""
A library for the antenna model
"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *

"""First import matlab data to python"""


class QP_Up(object):
    """
    This is for analyze the QP up transition rate with premeasurement as
    initialization.
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.Pre_to_RO_Time = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Single_Shot_Occupation_Filtered',
                             data_type2='RO_Pre_Drive_to_RO'):

        f_l = file_path[0]
        occ_2D = np.array([])
        self.Pre_to_RO_Time = np.array(noiselib.loadmat(f_l)[data_type2])
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))

        # print(occ_2D)
        # print(occ_2D.mean(axis=0))
        """Update parameters"""
        self.occ_1D_avg = occ_2D.mean(axis=0)

    def plot(self):
        time = self.Pre_to_RO_Time
        occ = self.occ_1D_avg

        plt.plot(time, occ, label='None')
        plt.xlabel('time (ns)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()


QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

# Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\P1Pre2021May14\Q4_P1_425mK\MATLABData

date = 'P1Pre2021May14'
experiment_name = ('Q4_P1_425mK')
# experiment_name = ('Q4_P1_375mK')
file_Number = np.arange(0, 30, 1)
# file_Number = np.arange(0, 12, 1)
QP_file = [QP_path.format(date, experiment_name, experiment_name) + '_{:03d}.mat'.format(i) for i in file_Number]
QP = QP_Up()
QP.add_data_from_matlab(QP_file)
QP.plot()
#
# QP_file = [QP_path]
# QP = QP_Recover()
# QP.add_data_from_matlab(QP_file)
# QP.plot()

def getGamma_pa(T, dfn=2e9, f0=120e9):
    """
    To calculate the theoretic photon assisted QP poisoning events based on
    the blackbody temperature and transmon's characteristic mode frequency
    :param T: temperature of the blackbody
    :param dfn: noise bandwidth
    :param f0: transmon's characteristic antenna mode frequency
    :return: Gamma_pa, the photon-assisted QP rate
    """
    Gamma_pa = dfn/(np.exp(h*f0/(k*T))-1)
    return Gamma_pa