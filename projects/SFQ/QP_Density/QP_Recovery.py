"""For analyze QP recover rate directly as Wang Chen's 2014 nature comminucation paper"""

import noiselib
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
from scipy.stats import linregress

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


class Phaseslip(object):

    def __init__(self):
        self.occ_2D = []
        self.poison_t = []
        self.idle_t = []
        self.phaseslip = []
        self.T1 = []
        self.x_QP = []
        self.n_QP = []
        self.C = self.getC()
        self.f_SFQ=2e9
        self.T1_R=40e-6
        self.T1_QP=100e-6

    def add_data_from_matlab(self, file_path,
                             data_type1='Weighted_Occupation',
                             data_type2='QB_Idle_Gate_Time',
                             data_type3='SFQ_Pulse_Duration'):

        for f in file_path:
            data = noiselib.loadmat(f)
            occ_2D = np.array(data[data_type1])
            idle_t = np.array(data[data_type2])
            poison_t = np.array(data[data_type3])

        """Update parameters and units"""
        # print('2D=', occ_2D)
        # print('time=', poison_t)
        self.occ_2D = occ_2D
        self.idle_t = idle_t*1e-9
        self.poison_t = poison_t*1e-9
        self.phaseslip = poison_t*1e-9*self.f_SFQ*3
        self.getT1()

    def getT1(self):
        T1 = []
        x_QP = []
        n_QP = []
        C = self.C
        occ2D = self.occ_2D
        Idle_array = self.idle_t[:50]
        for P1Idle in occ2D:
            # print('P1Idle=', P1Idle)
            P1_array = P1Idle[:50]
            t1 = self.fitToT1(P1_array, Idle_array)
            # print('t1=', t1)
            # T1.append(2e-6/(0.93-P1Idle[0]))
            n_qp = (1/t1 - 1/self.T1_R)*self.T1_QP
            # print('n_QP', n_QP)
            T1.append(t1)
            x_QP.append(1/(t1*C))
            n_QP.append(n_qp)
        self.T1 = np.array(T1)
        self.x_QP = np.array(x_QP)
        self.n_QP = np.array(n_QP)

    def fitToT1(self, P1_array, Idle_array):
        # T1 = 1
        # print('P1_array=', P1_array)
        # print('Idle_array=', Idle_array)
        slope, intercept, r_value, p_value, std_err = linregress(Idle_array, P1_array)
        # print('slope=', slope)
        T1 = -1/slope
        return T1

    def getC(self):
        Delta = 180*10**(-6)*e
        # print('Delta=', Delta)
        fq = 5*10**(9)
        # X = 8*fq*Delta/h
        # print('X=', X)
        C = np.sqrt(8*fq*Delta/h)
        # C = np.sqrt(X)
        # C = X**(1/2)
        # print('C=', C)
        return C

    def plot(self, mode='PhaseSlip_n_QP'):
        # Get common units
        poison_t  = self.poison_t * 10**6 # us
        # T1 = self.T1  # us
        T1 = self.T1 * 10**6 # us
        phase_slips = self.phaseslip
        x_QP = self.x_QP
        n_QP = self.n_QP

        if mode =='PoisonTime_T1':
            plt.plot(poison_t, T1, label='None')
            plt.xlabel('Poison time (us)')
            plt.ylabel('T1(us)')

        elif mode =='PhaseSlip_T1':
            plt.plot(phase_slips, T1, label='None')
            plt.xlabel('Phase Slips')
            plt.ylabel('T1(us)')

        elif mode =='x_QP_T1':
            plt.plot(phase_slips, x_QP, label='None')
            plt.xlabel('Phase Slips')
            plt.ylabel('x_QP')

        elif mode =='PhaseSlip_n_QP':
            plt.plot(phase_slips, n_QP, label='None')
            plt.xlabel('Phase Slips')
            plt.ylabel('n_QP')

        plt.grid()
        plt.legend()
        plt.show()


# QP_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM08/05-06-21'
#            '/T1_SFQ_Poison_Recover_Sweep/MATLABData'
#            '/T1_SFQ_Poison_Recover_Sweep_003_Weighted_Occupation_expfit.mat')
#
# QP_file = [QP_path]
# QP = QP_Recover()
# QP.add_data_from_matlab(QP_file)
# QP.plot()

# Z:\mcdermott-group\data\sfq\MCM_NIST\LIU\MCM08\05-06-21\T1_SFQ_Poison_Time_Sweep\MATLABData

phase_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM08/06-04-21'
           '/T1_SFQ_Poison_Time_Sweep/MATLABData'
           '/T1_SFQ_Poison_Time_Sweep_003.mat')

# Z:\mcdermott-group\data\sfq\MCM_NIST\LIU\MCM08\06-04-21\T1_SFQ_Poison_Time_Sweep\MATLABData

phase_file = [phase_path]
phase = Phaseslip()
phase.add_data_from_matlab(phase_file)
phase.plot()

