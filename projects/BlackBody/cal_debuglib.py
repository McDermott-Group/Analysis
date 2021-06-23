"""
A library for the antenna model
"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.constants import *
from scipy.stats import linregress
from scipy.stats import sem

from scipy.optimize import curve_fit


"""First import matlab data to python"""

class blobCenterStd(object):
    """
    Our P1 fluctuates with time. Is it from the measurement calibration
    Let's plot the IQ centers and std to see they fluctuate or not
    """

    def __init__(self):
        self.temp = None  # units mK
        self.gI = []
        self.gQ = []
        self.gIstd = []
        self.gQstd = []
        self.eI = []
        self.eQ = []
        self.eIstd = []
        self.eQstd = []

    def add_data_from_matlab(self, file_path, temp):
        gI = []
        gQ = []
        gIstd = []
        gQstd = []
        eI = []
        eQ = []
        eIstd = []
        eQstd = []
        for f in file_path:
            gI.append(noiselib.loadmat_ExptVars(f)['Ground_State_I'])
            gQ.append(noiselib.loadmat_ExptVars(f)['Ground_State_Q'])
            gIstd.append(noiselib.loadmat_ExptVars(f)['Ground_State_Istd'])
            gQstd.append(noiselib.loadmat_ExptVars(f)['Ground_State_Qstd'])
            eI.append(noiselib.loadmat_ExptVars(f)['Excited_State_I'])
            eQ.append(noiselib.loadmat_ExptVars(f)['Excited_State_Q'])
            eIstd.append(noiselib.loadmat_ExptVars(f)['Excited_State_Istd'])
            eQstd.append(noiselib.loadmat_ExptVars(f)['Excited_State_Qstd'])

        """Update parameters"""
        self.gI = np.array(gI)
        self.gQ = np.array(gQ)
        self.gIstd = np.array(gIstd)
        self.gQstd = np.array(gQstd)
        self.eI = np.array(eI)
        self.eQ = np.array(eQ)
        self.eIstd = np.array(eIstd)
        self.eQstd = np.array(eQstd)
        self.temp = temp

class TwoDFit_debug(object):
    """
    check how good is my earlier two dimensional fit. this is for single file
    """

    def __init__(self):
        self.Pe = 0.0
        self.params = []
        self.cov = []

        self.Is = []
        self.I = [] # avg Is
        self.Qs = []
        self.Q = [] # avg Qs

        self.gI = []
        self.gQ = []
        self.gIstd = []
        self.gQstd = []
        self.eI = []
        self.eQ = []
        self.eIstd = []
        self.eQstd = []

        self.theta = None

        self.IsQsProject2D = [] # on that line
        self.IsQsProject1D = [] # reduce the coordinated to 1d line
        self.Is_rot = []
        self.Qs_rot = []

        self.gI_rot= []
        self.gQ_rot = []
        self.eI_rot = []
        self.eQ_rot= []

    def add_data_from_matlab(self, file_path, temp,
                             data_type1='Is', data_type2='Qs',
                             data_type3='I', data_type4='Q'):
        f_l = file_path[0]
        Is = np.array(noiselib.loadmat(f_l)[data_type1])
        Qs = np.array(noiselib.loadmat(f_l)[data_type2])
        I = np.array(noiselib.loadmat(f_l)[data_type3]) # avg
        Q = np.array(noiselib.loadmat(f_l)[data_type4]) # avg

        self.gI.append(noiselib.loadmat_ExptVars(f_l)['Ground_State_I'])
        self.gQ.append(noiselib.loadmat_ExptVars(f_l)['Ground_State_Q'])
        self.gIstd.append(noiselib.loadmat_ExptVars(f_l)['Ground_State_Istd'])
        self.gQstd.append(noiselib.loadmat_ExptVars(f_l)['Ground_State_Qstd'])
        self.eI.append(noiselib.loadmat_ExptVars(f_l)['Excited_State_I'])
        self.eQ.append(noiselib.loadmat_ExptVars(f_l)['Excited_State_Q'])
        self.eIstd.append(noiselib.loadmat_ExptVars(f_l)['Excited_State_Istd'])
        self.eQstd.append(noiselib.loadmat_ExptVars(f_l)['Excited_State_Qstd'])

        """Update parameters and units"""
        self.Is = Is
        self.I = I
        self.Qs = Qs
        self.Q = Q
        self._getProjectionLine()
        self._getRotatedIQs()

    def _getProjectionLine(self):
        """
        I should just use coordinates rotation
        From the two centers, get the slope and interception of the line
        :return: the slope and interception
        """
        x1, y1 = self.gI[0], self.gQ[0]
        x2, y2 = self.eI[0], self.eQ[0]

        slope = (y2-y1)/(x2-x1)
        theta = np.arctan(slope)
        # interception = (x2*y1-x1*y2)/(x2-x1)
        # self.projection_line =[slope, interception]
        self.theta = theta

    def _getRotatedIQs(self):

        theta = self.theta
        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])

        Is_rot = []
        Qs_rot = []

        gI = self.gI
        gQ = self.gQ
        eI = self.eI
        eQ = self.eQ

        Is = np.array(self.Is).flatten()
        Qs = np.array(self.Qs).flatten()
        I = np.mean(self.I)
        Q = np.mean(self.Q)
        for i in range(len(Is)):
            I_rot, Q_rot= np.dot(np.array([Is[i]-I, Qs[i]-Q]), rot_matrix)
            Is_rot.append(I_rot)
            Qs_rot.append(Q_rot)

        gI_rot, gQ_rot = np.dot(np.array([gI[0]-I, gQ[0]-Q]), rot_matrix)
        eI_rot, eQ_rot = np.dot(np.array([eI[0]-I, eQ[0]-Q]), rot_matrix)


        self.gI_rot = [gI_rot]
        self.gQ_rot = [gQ_rot]
        self.eI_rot = [eI_rot]
        self.eQ_rot = [eQ_rot]

        self.Is_rot = np.array(Is_rot)
        self.Qs_rot = np.array(Qs_rot)

    def fit1D_Hist(self):
        """
        Fit the projected 1D hist to extract the occupation
        :return: The fitted Gaussian parameters
        """
        # def bimodal_Liu(x, A1, A2):
        #     return gauss(x, gI, gstd, A1) + gauss(x, eI, estd, A2)

        # Import the Hist and the centers

        gI = self.gI_rot[0]
        gstd = (self.gIstd[0]+self.gQstd[0])/2
        eI = self.eI_rot[0]
        estd = (self.eIstd[0]+self.eQstd[0])/2

        Is = self.Is_rot
        y, x, _ = plt.hist(Is, 100, density=1, alpha=0.5)
        x = (x[1:]+x[:-1])/2
        g_index = (np.abs(x-gI)).argmin()
        e_index = (np.abs(x-eI)).argmin()
        # print('y=', y)
        # print('x=', x)
        # print('g_index=', g_index)
        # print('e_index=', e_index)
        g_amp = np.mean(y[g_index-2: g_index+3])
        e_amp = np.mean(y[e_index-2: e_index+3])
        expected = (gI, gstd, g_amp, eI, estd, e_amp)
        # bounds = [(gI-gstd/10, gstd*0.8, g_amp*0.8, eI-estd/10, estd*0.8, e_amp*0.8),
        #           (gI+gstd/10, gstd*1.2, min(y, g_amp*1.2), eI+estd/10, estd*1.2, min(y, e_amp*1.2))]
        bounds = [(gI-gstd/10, gstd*0.8, g_amp*0.8, eI-estd/10, estd*0.8, e_amp*0.8),
                  (gI+gstd/10, gstd*1.2, min(max(y), g_amp*1.2), eI+estd/10, estd*1.2, min(max(y), e_amp*1.2))]
        params, cov = curve_fit(bimodal, x, y, bounds=bounds, p0=expected)

        # print('params=', params)
        gCenter, gSig, gAmp, eCenter, eSig, eAmp = \
            params[0], abs(params[1]), params[2], params[3], abs(params[4]), params[5]
        Pe = (eSig*eAmp)/(gSig*gAmp + eSig*eAmp)
        # print('Pe=', Pe)
        # print('Pe=', Pe)
        self.Pe = Pe
        self.params = params
        self.cov = cov

        # ##Show how good the fit is
        # plt.plot(x, gauss(x, gCenter, gSig, gAmp), label='g_fit')
        # plt.plot(x, gauss(x, eCenter, eSig, eAmp), label='e_fit')
        #
        # plt.plot(x, bimodal(x, *params), color='red', label='sum')
        # plt.legend()
        # plt.show()

    def plot_IQs(self):
        """
        plot IQ points
        :return:
        """
        Is = self.Is
        Qs = self.Qs
        plt.figure(1)
        plt.scatter(Is, Qs)
        # plt.show()

    def plot_rot_IQs(self):
        """
        plot IQ points
        :return:
        """
        Is = self.Is_rot
        Qs = self.Qs_rot
        plt.figure(2)
        plt.scatter(Is, Qs)
        # plt.show()

    def plot_rot_IQs_1D_Hist(self):
        """
        plot IQ points
        :return:
        """
        Is = self.Is_rot
        plt.figure(3)
        plt.hist(Is, bins=100)
        plt.show()

def gauss(x, mu, sigma, A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)
