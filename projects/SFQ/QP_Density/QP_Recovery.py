"""For analyze QP recover rate directly as Wang Chen's 2014 nature comminucation paper"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *

"""First import matlab data to python"""


class QP_Recover(object):

    def __init__(self):
        self.time_constant = []
        self.Gamma = []
        self.recover_time = []
        self.x_qp = []
        self.n_qp = []
        self.C = self.getC()

    def add_data_from_matlab(self, file_path,
                             data_type1='Extracted_Time_Constant',
                             data_type2='SFQ_Drive_to_QB'):
        for f in file_path:
            data = noiselib.loadmat_Liu(f)
            time_constant = 10 ** (-9) * np.array(data[data_type1])
            recover_time = 10 ** (-6) * np.array(data[data_type2])

        """Update parameters"""
        self.time_constant = time_constant
        self.recover_time = recover_time
        self.Gamma = 1/time_constant
        self.x_qp = self.Gamma/self.C   # with um^-3
        self.n_qp = self.Gamma/self.C

    def getC(self):
        Delta = 180*10**(-6)*e
        fq = 5*10**(9)
        C = (8*fq*Delta/h)**(1/2)
        return C

    def plot(self):
        Gamma = 1e-6 / self.time_constant
        t = 1e6 * self.recover_time
        x_qp = self.x_qp

        # plt.plot(t, x_qp, label='None')
        plt.plot(t, Gamma, label='None')
        plt.xlabel('time (us)')
        plt.ylabel('Gamma (1/us)')
        plt.yscale('log')
        # plt.ylabel('n_qp (1/m^3)')
        plt.grid()
        plt.legend()
        plt.show()


QP_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM08/05-06-21'
           '/T1_SFQ_Poison_Recover_Sweep/MATLABData'
           '/T1_SFQ_Poison_Recover_Sweep_003_Weighted_Occupation_expfit.mat')

QP_file = [QP_path]
QP = QP_Recover()
QP.add_data_from_matlab(QP_file)
QP.plot()
