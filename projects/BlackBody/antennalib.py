"""
A library for the antenna model
"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.stats import linregress

"""First import matlab data to python"""


class QP_Up(object):
    """
    This is for analyze the QP up transition rate with premeasurement as
    initialization.
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.Pre_to_RO_Time = []
        self.GammaUp = None
        self.temp = ''

    def add_data_from_matlab(self, file_path, temp,
                             data_type1='Weighted_Occupation_Filtered',
                             # data_type1='Single_Shot_Occupation_Filtered',
                             data_type2='RO_Pre_Drive_to_RO'):
        f_l = file_path[0]
        residual_reps = 0
        # residual_reps = noiselib.loadmat_ExptVars(f_l)['Residual_Reps']
        occ_2D = np.array([])
        self.Pre_to_RO_Time = np.array(
            noiselib.loadmat(f_l)[data_type2]) * 1e-9
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))

        for f in file_path[1:]:
            residual_reps += noiselib.loadmat_ExptVars(f)['Residual_Reps']
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            occ_2D = np.vstack((occ_2D, occ_1D))

        # print(occ_2D)
        # print(occ_2D.mean(axis=0))
        """Update parameters and units"""
        occ_1D_avg = np.array(occ_2D.mean(axis=0))
        error = [np.sqrt(p * (1 - p) / residual_reps) for p in occ_1D_avg]
        self.occ_1D_avg = occ_1D_avg
        self.occ_1D_avg_error_band = error
        self.fitToGammaUp()
        self.temp = temp

    def fitToGammaUp(self):
        time = self.Pre_to_RO_Time
        P1 = self.occ_1D_avg
        slope, intercept, r_value, p_value, std_err = linregress(time, P1)
        # print('slope=', slope)
        self.GammaUp = slope  # units in Hz
        self.occ_1D_avg_fit = slope * time + intercept
        return None

    def plot(self):
        time = self.Pre_to_RO_Time * 10 ** 6
        occ = self.occ_1D_avg
        occ_fit = self.occ_1D_avg_fit
        error = self.occ_1D_avg_error_band
        GammaUp = self.GammaUp

        plt.plot(time, occ_fit, 'o-', label=self.temp+'_GammaUp={0:.4g} Hz'.format(GammaUp))
        plt.plot(time, occ, 'k-', label=self.temp)
        plt.fill_between(time, occ - error, occ + error, alpha=0.2)
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.grid()
        plt.legend()
        plt.show()


def getGamma_pa(T, dfn=2e9, f0=120e9):
    """
    To calculate the theoretic photon assisted QP poisoning events based on
    the blackbody temperature and transmon's characteristic mode frequency
    :param T: temperature of the blackbody
    :param dfn: noise bandwidth
    :param f0: transmon's characteristic antenna mode frequency
    :return: Gamma_pa, the photon-assisted QP rate
    """
    Gamma_pa = dfn / (np.exp(h * f0 / (k * T)) - 1)
    return Gamma_pa
