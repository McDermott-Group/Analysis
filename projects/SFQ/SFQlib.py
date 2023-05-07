"""
A python library for SFQ project
Author: Chuan-Hong (Vincent) Liu
Date: 2021Dec28
"""
import noiselib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import rcParams
from scipy.constants import *
from scipy.stats import linregress
from scipy.stats import sem
from sklearn.mixture import GaussianMixture
from scipy.optimize import curve_fit
from scipy.linalg import expm
import random
from numpy import copy


def add_2Ddata_from_matlab(file_path, v1, v2, d):
    """

    :param file_path:
    :param v1: variable one
    :param v2: variable two
    :param d: dependent
    :return: v1, v2, and d arrays
    """
    data = noiselib.loadmat(file_path)
    datav1 = data[v1]
    datav2 = data[v2]
    datad = data[d]

    datad = np.rot90(datad)

    return datav1, datav2, datad

def add_CST_SimResultFromTxt(file):
    t = np.loadtxt(file, usecols=[0], skiprows=3)
    d = np.loadtxt(file, usecols=[1], skiprows=3)
    return [t, d]


class SFQ_3level_sweep(object):
    """
    From Robert's matlab code
    """

    def __init__(self):
        # qubit parameters
        self.w01 = 2*pi*4.903e9     # 01 transition
        self.w12 = 2*pi*4.625e9   # 12 transition
        self.c_qb = 77*1e-15    # qb capacitance
        self.c_sfq = 200*1e-18    # sfq-qb coupling capacitance
        self.t0 = 2 * (2*pi/self.w01)   # subharmonics * clock time

        # qubit states vectors
        self.g = np.array([1.0, 0, 0])  # ground state
        self.p = (1/np.sqrt(2))*np.array([1.0, 1.0, 0])  #
        self.p_i = (1/np.sqrt(2))*np.array([1.0, 1.0j, 0])  #
        self.m = (1/np.sqrt(2))*np.array([1.0, -1.0, 0])  #
        self.m_i = (1/np.sqrt(2))*np.array([1.0, -1.0j, 0])  #
        self.e = np.array([0.0, 1.0, 0])  # excited state

        self.proj1 = np.array([0.0, 1.0, 0])  # 1 state projection
        self.proj2 = np.array([0.0, 0.0, 1.0])  # 2 state projection

        #
        self.a = np.array([[0, 1.0, 0], [0, 0, np.sqrt(2)], [0, 0, 0]])
        self.a_dag = np.array([[0, 0, 0], [1.0, 0, 0], [0, np.sqrt(2), 0]])
        self.ufr = np.array([[1.0, 0, 0],
                             [0, np.exp(-1j*self.w01*self.t0), 0],
                             [0, 0, np.exp(-1j*(self.w01+self.w12)*self.t0)]])

    def analyze(self):
        # n = np.linspace(10, 500, 491)
        n = np.arange(10, 501, 1)
        # print(n)

        F_g = np.zeros(len(n))
        F_e = np.zeros(len(n))
        F_m = np.zeros(len(n))
        F_mi = np.zeros(len(n))
        F_p = np.zeros(len(n))
        F_pi = np.zeros(len(n))
        F_avg = np.zeros(len(n))

        p2_g = np.zeros(len(n))
        p2_e = np.zeros(len(n))
        p2_m = np.zeros(len(n))
        p2_mi = np.zeros(len(n))
        p2_p = np.zeros(len(n))
        p2_pi = np.zeros(len(n))
        p2_avg = np.zeros(len(n))

        ufr = 1.0*self.ufr
        proj2 = 1.0*self.proj2

        uni = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])  # unitary matrix

        for i in range(len(n)):

            psi_g = 1.0*self.g
            psi_e = 1.0*self.e
            psi_m = 1.0*self.m
            psi_mi = 1.0*self.m_i
            psi_p = 1.0*self.p
            psi_pi = 1.0*self.p_i

            phi = pi / (2.0 * n[i])
            # print('phi=', phi)
            usfq = expm((phi / 2) * (self.a - self.a_dag))
            # usfq = expm(self.a-self.a_dag)
            # print('usfq=', usfq)
            for ii in range(n[i]):
                psi_g = usfq.dot(psi_g)
                psi_e = usfq.dot(psi_e)
                psi_m = usfq.dot(psi_m)
                psi_mi = usfq.dot(psi_mi)
                psi_p = usfq.dot(psi_p)
                psi_pi = usfq.dot(psi_pi)
                uni = usfq.dot(uni)

                if ii < n[i]:   # free evolution
                    psi_g = ufr.dot(psi_g)
                    psi_e = ufr.dot(psi_e)
                    psi_m = ufr.dot(psi_m)
                    psi_mi = ufr.dot(psi_mi)
                    psi_p = ufr.dot(psi_p)
                    psi_pi = ufr.dot(psi_pi)
                    uni = ufr.dot(uni)

            F_g[i] = np.abs(np.dot(self.m, psi_g))**2
            F_e[i] = np.abs(np.dot(self.p, psi_e))**2
            F_m[i] = np.abs(np.dot(self.e, psi_m))**2
            F_mi[i] = np.abs(np.dot(self.m_i, psi_pi))**2
            F_p[i] = np.abs(np.dot(self.g, psi_p))**2
            F_pi[i] = np.abs(np.dot(self.p_i, psi_mi))**2
            F_avg[i] = (1/6.0)*(F_g[i]+F_e[i]+F_m[i]+F_mi[i]+F_p[i]+F_pi[i])

            p2_g[i] = np.abs(np.dot(proj2, psi_g))**2
            p2_e[i] = np.abs(np.dot(proj2, psi_e))**2
            p2_m[i] = np.abs(np.dot(proj2, psi_m))**2
            p2_mi[i] = np.abs(np.dot(proj2, psi_mi))**2
            p2_p[i] = np.abs(np.dot(proj2, psi_p))**2
            p2_pi[i] = np.abs(np.dot(proj2, psi_pi))**2
            p2_avg[i] = (1/6.0)*(p2_g[i]+p2_e[i]+p2_m[i]+p2_mi[i]+p2_p[i]+p2_pi[i])

        plt.plot(n, p2_g, 'b.', label='p2_g')
        plt.plot(n, p2_e, 'r-', label='p2_e')
        # plt.plot(n, p2_avg, 'r--')
        plt.plot(n, 1.0-F_avg, 'k--', label='avg error')
        # plt.plot(n, 1.0-F_g, label='avg error')
        # plt.plot(n, 1.0-F_e, label='avg error')
        # plt.plot(n, 1.0-F_m, label='F_m avg error')
        # plt.plot(n, 1.0-F_mi, label='F_mi avg error')
        # plt.plot(n, 1.0-F_p, label='F_p avg error')
        # plt.plot(n, 1.0-F_pi, label='F_pi avg error')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.show()


class SFQ_3level_sweep_error(object):
    """
    From Robert's matlab code
    """

    def __init__(self):
        # qubit parameters
        self.w01 = 2*pi*4.903e9     # 01 transition
        self.w12 = 2*pi*4.625e9   # 12 transition
        # self.w01 = 2*pi*5e9     # 01 transition
        self.w12 = 2*pi*4.8e9   # 12 transition
        self.c_qb = 77*1e-15    # qb capacitance
        self.c_sfq = 200*1e-18    # sfq-qb coupling capacitance
        self.t0 = 4 * (2*pi/self.w01)   # subharmonics * clock time

        # qubit states vectors
        self.g = np.array([1.0, 0, 0])  # ground state
        self.p = (1/np.sqrt(2))*np.array([1.0, 1.0, 0])  #
        self.p_i = (1/np.sqrt(2))*np.array([1.0, 1.0j, 0])  #
        self.m = (1/np.sqrt(2))*np.array([1.0, -1.0, 0])  #
        self.m_i = (1/np.sqrt(2))*np.array([1.0, -1.0j, 0])  #
        self.e = np.array([0.0, 1.0, 0])  # excited state

        self.proj0 = np.array([1.0, 0.0, 0])  # 0 state projection
        self.proj1 = np.array([0.0, 1.0, 0])  # 1 state projection
        self.proj2 = np.array([0.0, 0.0, 1.0])  # 2 state projection

        self.a = np.array([[0, 1.0, 0], [0, 0, np.sqrt(2)], [0, 0, 0]])
        self.a_dag = np.array([[0, 0, 0], [1.0, 0, 0], [0, np.sqrt(2), 0]])
        self.ufr = np.array([[1.0, 0, 0],
                             [0, np.exp(-1j*self.w01*self.t0), 0],
                             [0, 0, np.exp(-1j*(self.w01+self.w12)*self.t0)]])

    def analyze_3D(self):
        data_2D = []
        p2data_2D = []
        # p_SFQ = np.linspace(0.75, 1.0, 26)
        # p_SFQ = np.linspace(0.0, 0.05, 26)
        p_SFQ = np.linspace(0.96, 0.98, 6)
        for i in range(len(p_SFQ)):
            print('analyze_3D', p_SFQ[i])
            data = self._analyze_2D(P_SFQ=p_SFQ[i])
            data_1D = [p_SFQ[i], data[0][0], data[1][0]]
            p2data_1D = [p_SFQ[i], data[2][0], data[3][0]]
            data_2D.append(data_1D)
            p2data_2D.append(p2data_1D)
            # print('data_1D=', data_1D)
        # print('data_2D=', data_2D)
        # np.savetxt('SFQ_Error_LosePulses.txt', data_2D)
        # np.savetxt('SFQ_Error_GainPulses.txt', data_2D)
        # np.savetxt('SFQ_Error_GainPulses.txt', data_2D)
        # np.savetxt('SFQ_P2_GainPulses.txt', p2data_2D)
        # np.savetxt('SFQ_Error_LosePulses.txt', data_2D)
        # np.savetxt('SFQ_P2_LosePulses.txt', p2data_2D)
        data_2D = np.array(data_2D)
        p2data_2D = np.array(p2data_2D)
        plt.errorbar(p2data_2D[:, 0], p2data_2D[:, 1], yerr=p2data_2D[:, 2])
        plt.errorbar(data_2D[:, 0], data_2D[:, 1], yerr=data_2D[:, 2])
        plt.show()

    def _analyze_2D(self, P_SFQ = 1.0):
        n_SFQ_len = 75
        F_avg_2D_len = 1 # number of average
        n_SFQ = np.arange(0, n_SFQ_len, 1)
        F_avg_2D = np.zeros((F_avg_2D_len, n_SFQ_len))
        p2_avg_2D = np.zeros((F_avg_2D_len, n_SFQ_len))
        p_SFQ = P_SFQ
        # print(len(F_avg_2D))
        for i in range(len(F_avg_2D)):
            f_avg_1D = self._analyze_1D(n_SFQ_Len=n_SFQ_len, P_SFQ=p_SFQ)[0]
            p2_avg_1D = self._analyze_1D(n_SFQ_Len=n_SFQ_len, P_SFQ=p_SFQ)[1]
            F_avg_2D[i] = f_avg_1D
            p2_avg_2D[i] = p2_avg_1D

        F_avg_2D_avg = np.average(F_avg_2D, axis=0)
        F_avg_2D_se = np.std(F_avg_2D, axis=0)
        p2_avg_2D_avg = np.average(p2_avg_2D, axis=0)
        p2_avg_2D_se = np.std(p2_avg_2D, axis=0)
        # print('F_avg_2D_se=', F_avg_2D_se)
        max_ind = np.where(F_avg_2D_avg == F_avg_2D_avg.max())

        # print('max_ind=', max_ind)
        # print('1-F_avg_2D_avg[max_ind]=', (1-F_avg_2D_avg[max_ind])*1e3)
        # print('F_avg_2D_se[max_ind]=', 1e3*F_avg_2D_se[max_ind]/np.sqrt(F_avg_2D_len))
        # plt.plot(n_SFQ, 1-F_avg_2D_avg)
        # # plt.plot(n_SFQ, F_avg_2D_se)
        # # plt.xscale('log')
        # plt.yscale('log')
        # plt.legend()
        # plt.show()
        return 1-F_avg_2D_avg[max_ind], F_avg_2D_se[max_ind]/np.sqrt(F_avg_2D_len), p2_avg_2D_avg[max_ind]\
            , p2_avg_2D_se[max_ind]/np.sqrt(F_avg_2D_len)

    def _analyze_1D(self, n_SFQ_Len=100, P_SFQ=1.0):
        n_SFQ = np.arange(0, n_SFQ_Len, 1)   # number of SFQ pulses applied

        F_g = np.zeros(n_SFQ_Len)
        F_e = np.zeros(n_SFQ_Len)
        F_m = np.zeros(n_SFQ_Len)
        F_mi = np.zeros(n_SFQ_Len)
        F_p = np.zeros(n_SFQ_Len)
        F_pi = np.zeros(n_SFQ_Len)
        F_avg = np.zeros(n_SFQ_Len)

        p0_g = np.zeros(n_SFQ_Len)
        p0_e = np.zeros(n_SFQ_Len)
        p1_g = np.zeros(n_SFQ_Len)
        p1_e = np.zeros(n_SFQ_Len)
        p2_g = np.zeros(n_SFQ_Len)
        p2_e = np.zeros(n_SFQ_Len)

        p2_m = np.zeros(n_SFQ_Len)
        p2_mi = np.zeros(n_SFQ_Len)
        p2_p = np.zeros(n_SFQ_Len)
        p2_pi = np.zeros(n_SFQ_Len)
        p2_avg = np.zeros(n_SFQ_Len)

        ufr = 1.0*self.ufr
        proj0 = 1.0*self.proj0
        proj1 = 1.0*self.proj1
        proj2 = 1.0*self.proj2

        uni = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])  # unitary matrix
        n_pi = 96  # need 100 pulses to achieve pi gate
        phi = pi / n_pi
        usfq = expm((phi / 2) * (self.a - self.a_dag))

        psi_g = 1.0 * self.g
        psi_e = 1.0 * self.e
        psi_m = 1.0 * self.m
        psi_mi = 1.0 * self.m_i
        psi_p = 1.0 * self.p
        psi_pi = 1.0 * self.p_i

        p_SFQ = P_SFQ

        for i in n_SFQ:
            # print('i=', i)
            run_SFQ = np.random.choice(np.array([0, 1]), p=[1-p_SFQ, p_SFQ])
            # print('run_SFQ=', run_SFQ)
            if run_SFQ: # if run_SFQ==1, then apply the SFQ pulses, otherwise lose the pulse
                # H_SFQ
                psi_g = usfq.dot(psi_g)
                psi_e = usfq.dot(psi_e)
                psi_m = usfq.dot(psi_m)
                psi_mi = usfq.dot(psi_mi)
                psi_p = usfq.dot(psi_p)
                psi_pi = usfq.dot(psi_pi)
                uni = usfq.dot(uni)

            ## H_SFQ
            # psi_g = usfq.dot(psi_g)
            # psi_e = usfq.dot(psi_e)
            # psi_m = usfq.dot(psi_m)
            # psi_mi = usfq.dot(psi_mi)
            # psi_p = usfq.dot(psi_p)
            # psi_pi = usfq.dot(psi_pi)
            # uni = usfq.dot(uni)

            # H_fr
            psi_g = ufr.dot(psi_g)
            psi_e = ufr.dot(psi_e)
            psi_m = ufr.dot(psi_m)
            psi_mi = ufr.dot(psi_mi)
            psi_p = ufr.dot(psi_p)
            psi_pi = ufr.dot(psi_pi)
            uni = ufr.dot(uni)

            F_g[i] = np.abs(np.dot(self.m, psi_g))**2
            F_e[i] = np.abs(np.dot(self.p, psi_e))**2
            F_m[i] = np.abs(np.dot(self.e, psi_m))**2
            F_mi[i] = np.abs(np.dot(self.m_i, psi_pi))**2
            F_p[i] = np.abs(np.dot(self.g, psi_p))**2
            F_pi[i] = np.abs(np.dot(self.p_i, psi_mi))**2
            F_avg[i] = (1/6.0)*(F_g[i]+F_e[i]+F_m[i]+F_mi[i]+F_p[i]+F_pi[i])

            p0_g[i] = np.abs(np.dot(proj0, psi_g))**2
            p0_e[i] = np.abs(np.dot(proj0, psi_e))**2

            p1_g[i] = np.abs(np.dot(proj1, psi_g))**2
            p1_e[i] = np.abs(np.dot(proj1, psi_e))**2

            p2_g[i] = np.abs(np.dot(proj2, psi_g))**2
            p2_e[i] = np.abs(np.dot(proj2, psi_e))**2

            p2_m[i] = np.abs(np.dot(proj2, psi_m))**2
            p2_mi[i] = np.abs(np.dot(proj2, psi_mi))**2
            p2_p[i] = np.abs(np.dot(proj2, psi_p))**2
            p2_pi[i] = np.abs(np.dot(proj2, psi_pi))**2
            p2_avg[i] = (1/6.0)*(p2_g[i]+p2_e[i]+p2_m[i]+p2_mi[i]+p2_p[i]+p2_pi[i])

        # plt.plot(n_SFQ, p0_g+p1_g+p2_g, 'b.', label='p_g')
        # plt.plot(n_SFQ, p0_e+p1_e+p2_e, 'r-', label='p_e')

        # plt.plot(n_SFQ, p0_g, 'b.', label='p1_g')
        # plt.plot(n_SFQ, p0_e, 'r-', label='p1_e')

        # plt.plot(n_SFQ, p1_g, 'b.', label='p1_g')
        # plt.plot(n_SFQ, p1_e, 'r-', label='p1_e')
        #
        # plt.plot(n_SFQ, p2_g, 'b.', label='p2_g')
        # plt.plot(n_SFQ, p2_e, 'r-', label='p2_e')
        plt.plot(n_SFQ, 1.0-F_avg, 'k--', label='avg error')
        # plt.plot(n_SFQ, 1.0-F_g, label='avg error')
        # plt.plot(n_SFQ, 1.0-F_e, label='avg error')
        # plt.plot(n_SFQ, 1.0-F_m, label='F_m avg error')
        # plt.plot(n_SFQ, 1.0-F_mi, label='F_mi avg error')
        # plt.plot(n_SFQ, 1.0-F_p, label='F_p avg error')
        # plt.plot(n_SFQ, 1.0-F_pi, label='F_pi avg error')
        # plt.xscale('log')

        # plt.yscale('log')
        plt.legend()
        plt.show()
        return F_avg, p2_avg
        # return p2_g


class T1_QP_1D(object):
    """
    This is for extract 1D T1 value from the matlab data
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.QB_Idle_Gate_Time = []
        self.T1 = None
        self.Gamma_fit_parameters = []  # save each file's fitting parameteres
        self.fit_number_list = []
        self.T1_list = []
        self.params = None
        self.err = None
        self.fit_type = 'Exponential'

    def add_data_from_matlab(self, file_path,
                             data_type1='Weighted_Occupation',
                             # data_type1='Projected_Occupation',
                             data_type2='QB_Idle_Gate_Time',
                             fit_type='QP'):
        f_l = file_path[0]
        occ_2D = np.array([])
        Gamma_fit_parameters = np.array([])
        self.QB_Idle_Gate_Time = np.array(
            noiselib.loadmat(f_l)[data_type2])  # unit is ns
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        occ_2D = np.hstack((occ_2D, occ_1D))  # for future average use
        self.fit_type = fit_type

        if fit_type == 'QP':
            fit_params = self.fitToQPDecay(occ_1D)[1]
        Gamma_fit_parameters = np.hstack((Gamma_fit_parameters, fit_params))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])

            if fit_type == 'QP':
                fit_params = self.fitToQPDecay(occ_1D)[1]

            Gamma_fit_parameters = np.hstack(
                (Gamma_fit_parameters, fit_params))

            occ_2D = np.vstack((occ_2D, occ_1D))

        """Update parameters and units"""

        self.occ_1D_avg = occ_1D
        self.Gamma_fit_parameters = Gamma_fit_parameters

    # def fitLinear(self, occ_1D, pts):
    #     time = self.QB_Idle_Gate_Time[pts]
    #     P1 = occ_1D[pts]
    #     params, covariance = curve_fit(self.f_lin, time, P1, p0=[0.9, 30e4])
    #     self.params = params
    #     return params

    # def fitToExponential(self, occ_1D):
    #     """
    #     Extract the time constant, amp, offset and their uncertainties
    #     :param occ_1D:
    #     :return:
    #     """
    #     time = self.QB_Idle_Gate_Time
    #     P1 = occ_1D
    #     p0 =
    #     bounds = [(0.5, 1e4, 0), (1.05, 1e6, 0.4)]
    #     # params, covariance = curve_fit(self.f_exp_off, time, P1,
    #     #                                bounds=bounds, p0=[0.9, 1e5, 0.1])
    #     # params, covariance = curve_fit(self.f_exp_off, time, P1,
    #     #                                p0=[0.9, 2e4, 0.01])
    #     params, covariance = curve_fit(self.f_exp_off, time, P1, p0=p0, bounds=bounds)
    #     self.params = params
    #     return params

    def fitToQPDecay(self, occ_1D):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = self.QB_Idle_Gate_Time
        P1 = occ_1D
        p0 = [0.1, 20e3, 20e3, 0.95, 0.1]  # initial guess
        bounds = [(0, 1e3, 1e3, 0.5, 0), (10, 1e6, 1e6, 1.05, 0.5)]  # bounds
        params, covariance = curve_fit(self.f_QP, time, P1, p0=p0,
                                       bounds=bounds)
        self.params = params
        self.err = np.sqrt(
            np.diag(covariance))  # one standard deviation errors
        # print('params, err=', params, self.err)
        return params

    # def f_lin(self, t, intercept, Gamma):
    #     return intercept - t * Gamma

    # def f_exp(self, t, amp, Gamma):
    #     return amp * np.exp(-t * Gamma) + 0.05
    #
    # def f_exp_off(self, t, amp, Gamma, off):
    #     return amp * np.exp(-t * Gamma) + off

    def f_QP(self, t, n_qp, T_qp, T_1r, amp, off):
        return amp * np.exp(n_qp * (np.exp(-t / T_qp) - 1) - t / T_1r) + off

    def plot_histogram(self):
        Gamma = self.Gamma_fit_parameters
        T1 = [10 ** 6 / g for g in Gamma]  # in us
        plt.hist(T1, bins=10)
        plt.show()

    def plot(self):
        time = self.QB_Idle_Gate_Time
        occ = self.occ_1D_avg
        plt.figure(1)
        plt.plot(time, occ, 'k-')

        [n_qp, T_qp, T_1r, amp, off] = self.params
        [n_qp_std, T_qp_std, T_1r_std, amp_std, off_std] = self.err
        plt.plot(time, self.f_QP(time, n_qp, T_qp, T_1r, amp, off))

        n_qp_str = 'n_qp={0:.4g} $\pm$ {1:.4g}'.format(n_qp, n_qp_std)
        T_qp_str = 'T_qp={0:.4g} $\pm$ {1:.4g} ns'.format(T_qp, T_qp_std)
        T_1r_str = 'T1_r={0:.4g} $\pm$ {1:.4g} ns'.format(T_1r, T_1r_std)

        plt.xlabel('time (ns)')
        plt.ylabel('P1')
        # plt.yscale('log')
        plt.title(n_qp_str + '\n' + T_qp_str + ', ' + T_1r_str)
        plt.grid()
        plt.legend()
        plt.show()


class T1_QP_Paper(object):
    def __init__(self):
        self.P1_2D = None
        self.params_2D = None
        self.std_2D = None

        self.QB_Idle_Gate_Time = None
        self.sweep_variable_list = None
        self.sweep_variable_name = None

        self.T1_1D_t = None
        self.T1_1D_noPoison = None
        self.T1_1D_noPoisonFit = None
        self.T1_1D_Poison = None
        self.T1_1D_PoisonFit = None

    def add_data_from_matlab(self, file_path,
                             data_type='Weighted_Occupation',
                             # data_type='Projected_Occupation',
                             data_Idle_Time='QB_Idle_Gate_Time',
                             data_sweep_variable='SFQ_Pulse_Duration'
                             ):
        f = file_path[0]
        data = noiselib.loadmat(f)
        P1_2D = np.array(data[data_type])
        QB_Idle_Gate_Time = np.array(data[data_Idle_Time])
        sweep_variable_list = np.array(data[data_sweep_variable])

        self.P1_2D = P1_2D

        self.QB_Idle_Gate_Time = QB_Idle_Gate_Time
        self.sweep_variable_list = sweep_variable_list
        self.sweep_variable_name = data_sweep_variable

        self._analyze()
        # self._plot()
        return

    def _analyze(self):
        P1_2D = self.P1_2D
        t = self.QB_Idle_Gate_Time
        params_2D = []
        std_2D = []
        ### 1D fit
        for i in range(len(P1_2D)):
            # print('i=', i)
            P1_1D = P1_2D[i]
            params, std = self._fitToQPDecay(P1_1D, t)
            params_2D.append(params)
            std_2D.append(std)
            if i == 0:
                print('i=', i)
                print('params=', params)
                # print('error=', np.sqrt(np.diag(std))[0])
                print('params=', "n_qp, T_qp, T_1r, amp, off")
                self.T1_1D_t = t
                self.T1_1D_noPoison = P1_1D
                self.T1_1D_noPoisonFit = self._f_QP(t, *params)

            # if i > 4 and i < 21:
            #     # self.T1_1D_Poison = P1_1D
            #     # self.T1_1D_PoisonFit = self._f_QP(t, *params)
            #     plt.scatter(t, self.T1_1D_noPoison, label=str(i))
            #     plt.plot(t, self.T1_1D_noPoisonFit)
            #     plt.scatter(t, P1_1D)
            #     plt.plot(t, self._f_QP(t, *params))
            #     plt.legend()
            #     plt.show()

            if i == 7:
                print('i=', i)
                print('params=', params)
                # print('error=', np.sqrt(np.diag(std))[0])
                print('params=', "n_qp, T_qp, T_1r, amp, off")
                self.T1_1D_Poison = P1_1D
                self.T1_1D_PoisonFit = self._f_QP(t, *params)

            # if i == 20:
            #     print('i=', i)
            #     print('params=', params)
            #     print('params=', "n_qp, T_qp, T_1r, amp, off")
            #     plt.plot(t, P1_1D)
            #     plt.plot(t, self._f_QP(t, *params))
            #     plt.show()

        ### 2D data update
        self.params_2D = np.array(params_2D)
        self.std_2D = np.array(std_2D)
        return

    def _plot(self):
        params_2D = self.params_2D
        std_2D = self.std_2D
        sweep_variable_list = self.sweep_variable_list
        sweep_variable_name = self.sweep_variable_name
        # print('params_2D=', params_2D)
        numberofJJs = 4
        CW_freq = 1.200  # units in GHz
        # ps = numberofJJs * CW_freq  # phase slips conversion, 3 is the junction number,
        ps = 1
        plt.figure()
        plt.errorbar(sweep_variable_list * ps, params_2D[:, 0],
                     yerr=std_2D[:, 0], ls='none')
        plt.scatter(sweep_variable_list * ps, params_2D[:, 0])

        print(repr(params_2D[:, 0]))
        plt.show()

        return

    def _fitToQPDecay(self, P1_1D, time):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = time
        P1 = P1_1D
        p0 = [0.2, 6e3, 20e3, 0.95, 0.1]  # initial guess
        bounds = [(0, 5e3, 15e3, 0.5, 0), (2, 9e3, 30e3, 1.05, 0.3)]  # bounds
        params, covariance = curve_fit(self._f_QP, time, P1, p0=p0,
                                       bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        return list(params), list(std)

    def _f_QP(self, t, n_qp, T_qp, T_1r, amp, off):
        return amp * np.exp(n_qp * (np.exp(-t / T_qp) - 1) - t / T_1r) + off


class T1_QP_2D(object):
    def __init__(self):
        self.P1_2D = None
        self.params_2D = None
        self.std_2D = None

        self.QB_Idle_Gate_Time = None
        self.sweep_variable_list = None
        self.sweep_variable_name = None

    def add_data_from_matlab(self, file_path,
                             data_type='Weighted_Occupation',
                             # data_type='Projected_Occupation',
                             data_Idle_Time='QB_Idle_Gate_Time',
                             data_sweep_variable='SFQ_Pulse_Duration'
                             # data_sweep_variable='SFQ_Drive_to_QB'
                             ):
        f = file_path[0]
        data = noiselib.loadmat(f)
        P1_2D = np.array(data[data_type])
        QB_Idle_Gate_Time = np.array(data[data_Idle_Time])
        sweep_variable_list = np.array(data[data_sweep_variable])

        self.P1_2D = P1_2D

        self.QB_Idle_Gate_Time = QB_Idle_Gate_Time
        self.sweep_variable_list = sweep_variable_list
        self.sweep_variable_name = data_sweep_variable

        self._analyze()
        # self._plot()
        return

    def _analyze(self):
        P1_2D = self.P1_2D
        t = self.QB_Idle_Gate_Time
        params_2D = []
        std_2D = []
        ### 1D fit
        for i in range(len(P1_2D)):
            P1_1D = P1_2D[i]
            params, std = self._fitToQPDecay(P1_1D, t)
            params_2D.append(params)
            std_2D.append(std)
            # if i == 0:
            #     print('params=', params)
            #     print('params=', "n_qp, T_qp, T_1r, amp, off")
            #     plt.plot(t, P1_1D)
            #     plt.plot(t, self._f_QP(t, *params))
            #     plt.show()

        ### 2D data update
        self.params_2D = np.array(params_2D)
        self.std_2D = np.array(std_2D)
        return

    def _plot(self):
        params_2D = self.params_2D
        std_2D = self.std_2D
        sweep_variable_list = self.sweep_variable_list
        sweep_variable_name = self.sweep_variable_name
        # print('params_2D=', params_2D)
        numberofJJs = 4
        CW_freq = 1.200  # units in GHz
        ps = numberofJJs * CW_freq  # phase slips conversion, 3 is the junction number,
        # ps = 1
        plt.figure()
        plt.errorbar(sweep_variable_list * ps, params_2D[:, 0],
                     yerr=std_2D[:, 0], ls='none')
        plt.scatter(sweep_variable_list * ps, params_2D[:, 0])
        print(repr(params_2D[:, 0]))
        if 0:  # recovery rate fitting
            r = 22
            t = sweep_variable_list[:r]
            n_qp = params_2D[:, 0][:r]
            params, std = self._fitToExp(n_qp, t)
            [T, amp, off] = params
            plt.plot(t, self._f_exp_off(t, *params),
                     label='recover time={0:.4g} (ns)'.format(T))
            plt.xlabel('Recover time (ns)')

        if 1:  # n_qp fitting
            ps = sweep_variable_list * 9
            rate = 6e-5
            [m, b] = np.polyfit(ps, params_2D[:, 0], 1)
            # plt.plot(ps, ps*m + b, label='QPs/slip={0:.4g}'.format(m), c='k')
            plt.xlabel('SFQ phase slips', fontsize=14)
        # plt.tick_params(labelsize=14)
        # plt.ylabel('$n_{\mathrm{QP}}$', fontsize=14)
        # # plt.xticks([0, 4000, 8000, 12000, 16000, 20000])
        # # plt.ylim([0, 2])
        # plt.legend(frameon=False, loc=2, prop={'size': 14})
        # plt.show()
        return

    def _fitToQPDecay(self, P1_1D, time):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param occ_1D:
        :return:
        """
        time = time
        P1 = P1_1D
        p0 = [0.2, 6e3, 20e3, 0.95, 0.1]  # initial guess
        bounds = [(0, 5e3, 15e3, 0.5, 0), (2, 9e3, 30e3, 1.05, 0.3)]  # bounds
        params, covariance = curve_fit(self._f_QP, time, P1, p0=p0,
                                       bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        return list(params), list(std)

    def _f_QP(self, t, n_qp, T_qp, T_1r, amp, off):
        return amp * np.exp(n_qp * (np.exp(-t / T_qp) - 1) - t / T_1r) + off

    def _fitToExp(self, n_qp, time):
        """
        extract recovery rate
        :param n_qp:
        :param time:
        :return:
        """
        time = time
        n_qp = n_qp
        p0 = [20e3, 2, 0]  # initial guess
        bounds = [(1e3, 0, 0), (1e5, 10, 1)]  # bounds
        params, covariance = curve_fit(self._f_exp_off, time, n_qp, p0=p0,
                                       bounds=bounds)
        std = np.sqrt(np.diag(covariance))  # one standard deviation errors
        return list(params), list(std)

    def _f_exp_off(self, t, T, amp, off):
        return amp * np.exp(-t / T) + off


class T1_QP_2D_Linear(object):
    def __init__(self):
        self.P1_2D = None

        self.T1_2D_stats = None
        self.T1_1D_avg = None
        self.T1_1D_se = None
        self.NoF = 1    # number of files

        self.params_2D = None
        self.std_2D = None

        self.QB_Idle_Gate_Time = None
        self.sweep_variable_list = None
        self.sweep_variable_name = None

    def add_data_from_matlab(self, file_path,
                             data_P1_2D='Weighted_Occupation',
                             data_Idle_Time='QB_Idle_Gate_Time',
                             data_sweep_variable='SFQ_Pulse_Duration'
                             ):


        f0= file_path[0]
        data = noiselib.loadmat(f0)
        P1_2D = np.array(data[data_P1_2D])
        QB_Idle_Gate_Time = np.array(data[data_Idle_Time]) / 1e9
        sweep_variable_list = np.array(data[data_sweep_variable])
        self.QB_Idle_Gate_Time = QB_Idle_Gate_Time
        self.sweep_variable_list = sweep_variable_list
        self.sweep_variable_name = data_sweep_variable
        T1_1D = self._analyze(P1_2D)
        T1_2D_stats = T1_1D

        for f in file_path[1:]:
            self.NoF = self.NoF + 1
            data = noiselib.loadmat(f)
            P1_2D = np.array(data[data_P1_2D])
            T1_1D = self._analyze(P1_2D)
            T1_2D_stats = np.vstack((T1_2D_stats, T1_1D))

        self.T1_2D_stats = T1_2D_stats
        self._T1Stats()
        return None

    def _analyze(self, P1_2D):
        t = self.QB_Idle_Gate_Time
        params_2D = []
        std_2D = []
        ### 1D fit
        for i in range(len(P1_2D)):
            P1_1D = P1_2D[i]
            params, std = self._fitToLinear(P1_1D, t)
            params_2D.append(params)
            std_2D.append(std)

        ### 2D data update
        params_2D = np.asarray(params_2D)
        std_2D = np.asarray(std_2D)
        self.params_2D = params_2D
        self.std_2D = std_2D
        # print('params_2D[:, 0]=', params_2D[:, 0])
        return params_2D[:, 0]

    def _T1Stats(self):
        T1_2D_stats = self.T1_2D_stats
        NoF = self.NoF
        if self.NoF == 1:
            self.T1_1D_avg = T1_2D_stats
        else:   # more than one file is processed
            self.T1_1D_avg = np.average(self.T1_2D_stats, axis=0)
            T1_1D_std = np.std(self.T1_2D_stats, axis=0)
            self.T1_1D_se = T1_1D_std/np.sqrt(NoF)

        return None

    def plot(self):
        # params_2D = self.params_2D
        # std_2D = self.std_2D
        sweep_variable_list = self.sweep_variable_list
        sweep_variable_name = self.sweep_variable_name
        # numberofJJs = 4
        f_on = 1.22617  # on resonant drive freq
        # f_off = 1.200  # off resonant drive freq
        # f_off = 1.0810  # off resonant drive freq
        f_off = 1.21  # off resonant drive freq
        avgClen_on = 91   # ns avg clifford gate length = 91 ns on resonance
        avgClen_off = (f_on/f_off)*avgClen_on   # ns off resonance
        if self.NoF == 1:
            plt.plot(sweep_variable_list/avgClen_off, self.T1_1D_avg)
        else:
            plt.errorbar(sweep_variable_list/avgClen_off, self.T1_1D_avg, yerr=self.T1_1D_se)

        plt.xlabel('Equivalent Clifford sequence length')
        plt.ylabel('T1 (us)')
        plt.title('T1 vs Off resonant poisoning')
        plt.show()
        return

    def _fitToLinear(self, P1_1D, time):
        """
        extract recovery rate
        :param n_qp:
        :param time:
        :return:
        """
        n = 20  # first 20 elements
        time = time[:n]
        P1_1D = P1_1D[:n]
        # p0 = [20e3, 2, 0]  # initial guess
        # bounds = [(1e3, 0, 0), (1e5, 10, 1)]  # bounds
        [[gamma, b], cov] = np.polyfit(time, P1_1D, 1, cov=True)
        std = np.sqrt(np.diag(cov))  # one standard deviation errors
        return list([-1e6 / gamma, b]), list(std)  # converts to T1 with us as
        # units


class T2_1D(object):
    """
    This is for extract 1D T1 value from the matlab data
    """

    def __init__(self):
        self.occ_1D_avg = []
        self.occ_1D_avg_fit = []
        self.occ_1D_avg_error_band = []
        self.QB_Idle_Gate_Time = []
        self.T2 = None
        self.fit_parameters = []  # save each file's fitting parameteres
        self.fit_number_list = []
        self.T2_list = []
        self.params = None

    def add_data_from_matlab(self, file_path,
                             # data_type1='Weighted_Occupation',
                             data_type1='Projected_Occupation',
                             data_type2='QB_Idle_Gate_Time'):
        f_l = file_path[0]
        occ_2D = np.array([])
        fit_parameters = np.array([])
        self.QB_Idle_Gate_Time = np.array(
            noiselib.loadmat(f_l)[data_type2])
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        fit_params = self.fitToRamsey(occ_1D, f_l[-7:-4])[3]
        # print('f_l[-7:-4])=', f_l[-7:-4])
        fit_parameters = np.hstack((fit_parameters, fit_params))

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            fit_params = self.fitToRamsey(occ_1D, f[-7:-4])[3]
            # print('f[-7:-4])=', f[-7:-4])
            fit_parameters = np.hstack(
                (fit_parameters, fit_params))

        """Update parameters and units"""

        self.occ_1D_avg = occ_1D
        self.fit_parameters = np.array(fit_parameters)

    def fitToRamsey(self, occ_1D, file_number):
        """
        Extract the time constant, amp, offset and their uncertainties
        :param file_number:
        :param occ_1D:
        :return:
        """
        time = self.QB_Idle_Gate_Time
        P1 = occ_1D
        bounds = [(0.0, 1e4, 1e4, 500, 0, 0.0), (1.0, 1e6, 1e6, 2000, pi, 0.8)]
        p0 = [0.8, 10e4, 10e4, 1000, 0, 0.5]
        # params, covariance = curve_fit(self.f_Ramsey, time, P1, bounds=bounds,
        #                                p0=p0)
        params, covariance = curve_fit(self.f_Ramsey, time, P1,
                                       p0=p0)
        self.occ_1D_avg = occ_1D
        self.params = params
        self._plot_fit(file_number, save=False)

        return params

    def f_Ramsey(self, t, amp, T_white, T_corr, T_Ramsey, phase_offset,
                 amp_offset):
        return amp * np.exp(-t / T_white - (t / T_corr) * (t / T_corr)) * \
               np.sin(2 * pi * (t / T_Ramsey) + phase_offset) + amp_offset

    def _plot_fit(self, file_number, save=False):
        time = self.QB_Idle_Gate_Time
        occ = self.occ_1D_avg
        [amp, T_white, T_corr, T_Ramsey, phase_offset,
         amp_offset] = self.params
        plt.figure()
        plt.plot(time, occ, 'k-')
        plt.plot(time, [amp * np.exp(-t / T_white - (t / T_corr) ** 2) * \
                        np.sin(2 * pi * (
                                t / T_Ramsey) + phase_offset) + amp_offset
                        for t in time], '--', label='fit')
        plt.xlabel('time (ns)')
        plt.ylabel('P1')
        plt.grid()
        # amp_str = 'amp={0:.4g} '.format(amp)
        T_white_str = 'T_whilte={0:.4g} (ns) '.format(T_white)
        T_corr_str = 'T_corr={0:.4g} (ns) '.format(T_corr)
        T_Ramsey_str = 'T_Ramsey={0:.4g} (ns)'.format(T_Ramsey)
        plt.title(T_white_str + '\n' + T_corr_str + T_Ramsey_str)
        plt.legend()
        if save:
            plt.savefig(file_number + '.png')
        else:
            plt.show()
        # plt.close()

    def plot_Detuning(self):
        Detuning = 1e3 / (self.fit_parameters)  # units in MHz
        x_len = len(Detuning)
        # sfq_pulse = np.linspace(0, 1000, 51) * 3 * 3
        # sfq_pulse = np.linspace(0, 5000, 500) * 1.2355 * 3
        # sfq_pulse = np.linspace(0, 6500, 13)
        # sfq_pulse = np.linspace(0, 500, 2)
        unit = 1000
        sfq_pulse = np.linspace(0, unit * (x_len - 1), x_len)
        plt.figure()
        plt.plot(sfq_pulse, Detuning)
        # plt.xlabel('phase slips (clock time reference)')
        plt.xlabel('Pulse Time (ns)')
        plt.ylabel('Detuning (MHz)')
        # plt.title('2021Dec25_Drive6.3225GHz')
        tt = '2022May19_1.2355GHz_Q1'
        plt.title(tt)
        # plt.savefig(tt + '.png')
        plt.show()


class QP_Up(object):
    """
    This is for analyze the QP up transition rate
    """

    def __init__(self):
        self.occ_1D = []
        self.occ_1D_fit = []
        self.SFQ_freq = None  # units is GHz
        # self.occ_1D_error_band = []
        self.Time_Dep = []
        self.GammaUp = None
        self.GammaUp_std_err = None
        self.temp = None  # units mK
        self.GammaUp_list = []

    def add_data_from_matlab(self, file_path,
                             data_type1='Projected_Occupation',
                             data_type2='SFQ_Pulse_Duration'):
        f_l = file_path[0]
        occ_2D = np.array([])

        # SFQ units is GHz
        SFQ_freq = noiselib.loadmat_ExptVars(f_l)['SFQ_Drive_Frequency']
        # print('SFQ_freq=', SFQ_freq)

        Time_Dep = np.array(noiselib.loadmat(f_l)[data_type2]) * 1e-9
        occ_1D = np.array(noiselib.loadmat(f_l)[data_type1])
        Time_Dep = Time_Dep[:]
        occ_1D = occ_1D[:]
        occ_2D = np.hstack((occ_2D, occ_1D))
        #
        # for f in file_path[1:]:
        #     data = noiselib.loadmat(f)
        #     occ_1D = np.array(data[data_type1])
        #     occ_2D = np.vstack((occ_2D, occ_1D))
        """Update parameters and units"""

        self.Time_Dep = Time_Dep
        self.SFQ_freq = SFQ_freq
        self.occ_1D = occ_1D
        self.fitToGammaUp()

    def fitToGammaUp(self):
        time = self.Time_Dep
        P1 = self.occ_1D
        slope, intercept, r_value, p_value, std_err = linregress(time, P1)
        self.GammaUp = slope  # units in Hz
        self.occ_1D_fit = slope * time + intercept
        self.GammaUp_std_err = std_err
        return None

    def plot(self, name='', save=False):
        time = self.Time_Dep * 10 ** 6
        occ = self.occ_1D
        occ_fit = self.occ_1D_fit
        GammaUp = self.GammaUp
        SFQ_freq = self.SFQ_freq  # converted in to GHz
        phase_slips_rate = SFQ_freq * 3 * 1e9  # 3 junctions phase slips
        plt.figure()
        plt.plot(time, occ_fit, 'o-',
                 label='GammaUp={0:.4g}Hz'.format(GammaUp))
        plt.plot(time, occ, 'k-')
        plt.xlabel('time (us)')
        plt.ylabel('P1')
        # plt.yscale('log')

        # n_qp_str = 'n_qp={0:.4g} $\pm$ {1:.4g}'.format(n_qp, n_qp_std)
        # T_qp_str = 'T_qp={0:.4g} $\pm$ {1:.4g} ns'.format(T_qp, T_qp_std)
        # T_1r_str = 'T1_r={0:.4g} $\pm$ {1:.4g} ns'.format(T_1r, T_1r_std)
        SFQ_str = 'SFQ Drive Freq = {0:.4g} GHz'.format(SFQ_freq)
        Up_phase_slip_str = 'Induced P1 per phase slip = {0:.4g}'.format(
            GammaUp / phase_slips_rate)

        plt.title(SFQ_str + '\n' + Up_phase_slip_str)

        plt.grid()
        plt.legend()
        plt.show()
        return


class DeepSubharmonics(object):
    """
    This is for average multi SFQ-based rabi to clearly see the steps from
    discretization of SFQ operation
    """

    def __init__(self):
        self.occ_1D = []
        self.occ_1D_avg = []
        self.occ_1D_std = []
        self.Time_Dep = []

    def add_data_from_matlab(self, file_path,
                             # data_type1='I',
                             # data_type1='Q',
                             # data_type1='Amplitude',
                             # data_type1='Phase',
                             # data_type1='Projected_Occupation',
                             data_type1='Weighted_Occupation',
                             # data_type2='SFQ_Pulse_Duration'):
                             data_type2='SFQ_Drive_to_RO'):
        f0 = file_path[0]
        data0 = noiselib.loadmat(f0)
        occ_2D = data0[data_type1]
        Time_Dep = data0[data_type2]
        icheck = 0
        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            occ_1D = np.array(data[data_type1])
            # print('occ_1D[0]=', occ_1D[0])
            # if occ_1D[0] > 0.1 and occ_1D[0] < 0.22:
            occ_2D = np.vstack((occ_2D, occ_1D))
                # icheck = icheck + 1
        # print('icheck=', icheck)

        """Update parameters"""
        self.Time_Dep = Time_Dep
        self.occ_1D_avg = np.average(occ_2D, axis=0)
        self.occ_1D_std = np.std(occ_2D, axis=0)/np.sqrt(len(file_path))

    def plot(self, name='', save=False):
        time = self.Time_Dep
        occ = self.occ_1D_avg
        std = self.occ_1D_std
        plt.figure()
        # plt.plot(time, occ, 'k-')
        plt.errorbar(time, occ, yerr=std)
        plt.xlabel('time (ns)', fontsize=24)
        plt.ylabel('P1', fontsize=24)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        # SFQ_str = 'SFQ Drive Freq = {0:.4g} GHz'.format(SFQ_freq)
        # Up_phase_slip_str = 'Induced P1 per phase slip = {0:.4g}'.format(GammaUp/phase_slips_rate)

        # plt.title(SFQ_str + '\n' + Up_phase_slip_str)

        plt.grid()
        plt.legend()
        plt.show()
        return


class RB(object):
    """
    Analyze the qubit gate fidelity. Import population data, fit, and
    extract the gate fidelity.
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords
        self.NIP = np.array([])  # no interleaved
        self.IP = np.array([])  # interleaved

        self.NIP_avg = np.array([])  # no interleaved avg
        self.NIP_std = np.array([])  # no interleaved standard deviation
        self.NIP_ste = np.array([])  # no interleaved standard error
        self.IP_avg = np.array([])  # interleaved avg
        self.IP_std = np.array([])  # interleaved standard deviation
        self.IP_ste = np.array([])  # interleaved standard error

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             data_type1='No_Interleaved_Probability',
                             data_type2='Interleaved_Probability'):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        self.IntGate = noiselib.loadmat_ExptVars(f0)['Interleaving_Gate']
        NIP_2D = data0[data_type1]
        IP_2D = data0[data_type2]

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            NIP_1D = np.array(data[data_type1])
            IP_1D = np.array(data[data_type2])
            NIP_2D = np.vstack((NIP_2D, NIP_1D))
            IP_2D = np.vstack((IP_2D, IP_1D))

        self.NIP = NIP_2D
        self.IP = IP_2D

    def data_analysis(self):
        self._convert_data()
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _convert_data(self):
        """
        Make the Number of cliffords ascending order and also the
        probability to fidelity
        :return:
        """
        self.NoC = self.NoC
        m = len(self.NoC)
        l1 = len(self.NIP)
        for i in range(l1):
            self.NIP[i] = np.ones(m) - self.NIP[i]
            self.IP[i] = np.ones(m) - self.IP[i]

        self.NIP_avg = np.average(self.NIP, axis=0)
        self.NIP_std = np.std(self.NIP, axis=0)
        self.NIP_ste = self.NIP_std / np.sqrt(l1)
        self.IP_avg = np.average(self.IP, axis=0)
        self.IP_std = np.std(self.IP, axis=0)
        self.IP_ste = self.IP_std / np.sqrt(l1)

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        NIP_avg = self.NIP_avg
        IP_avg = self.IP_avg
        NoC = self.NoC
        p0 = [0.5, 0.985, 0.5]
        params_NIP, covariance_NIP = curve_fit(self._function, NoC, NIP_avg,
                                               p0=p0)
        params_IP, covariance_IP = curve_fit(self._function, NoC, IP_avg,
                                             p0=p0)
        self.params_NIP = params_NIP
        self.params_IP = params_IP
        self.err_NIP = np.sqrt(np.diag(covariance_NIP))
        self.err_IP = np.sqrt(np.diag(covariance_IP))
        # print('params, err=', params_NIP, self.err_NIP)
        # print('params, err=', params_IP, self.err_IP)

    def _function(self, NoC, A, p, B):
        return A * p ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the (non) interleaved seq error and fidelity, and
        find out the interleaved gate's fidelity
        :return:
        """
        p_NIP = self.params_NIP[1]
        p_IP = self.params_IP[1]
        err_NIP = self.err_NIP[1]
        err_IP = self.err_IP[1]

        r_NI = (1 - p_NIP) / 2  # non-interleaved sequence error
        r_I = (1 - p_IP) / 2  # interleaved sequence error

        r_int = (1 - p_IP / p_NIP) / 2.0  # interleaved gate error
        F = 1 - r_int  # interleaved gate fidelity
        F_std = F * np.sqrt((err_NIP / p_NIP) ** 2 + (err_IP / p_IP) ** 2)
        self.r_NI = r_NI
        self.r_I = r_I
        self.F = F
        self.F_std = F_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)
        IntGate = self.IntGate
        NIP_avg = self.NIP_avg
        NIP_ste = self.NIP_ste
        IP_avg = self.IP_avg
        IP_ste = self.IP_ste
        """fitting results"""
        params_NIP = self.params_NIP
        params_IP = self.params_IP
        r_NI = self.r_NI
        r_I = self.r_I
        F = self.F
        F_std = self.F_std

        """real plot"""
        plt.figure()
        plt.errorbar(
            NoC, NIP_avg, yerr=NIP_ste, fmt="o", color='k', capsize=3.0,
            ms=5.0,
            label='Interleaved Gate: {}, AvgError = {:.3f} %'.format('None',
                                                                     r_NI * 100))
        plt.plot(NoC_fit, self._function(NoC_fit, *params_NIP), 'k')

        plt.errorbar(
            NoC, IP_avg, yerr=IP_ste, fmt="o", color='r', capsize=3.0, ms=5.0,
            label='Interleaved Gate: {}, AvgError = {:.3f} %'.format(IntGate,
                                                                     r_I * 100))
        plt.plot(NoC_fit, self._function(NoC_fit, *params_IP), 'r')

        plt.title('{} Gate Fidelity = {:.3f} $\pm$ {:.3f} %'.
                  format(IntGate, F * 100, F_std * 100))

        plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Sequence Fidelity')
        plt.show()


class Purity_Paper(object):
    """
    Analyze the qubit gate incoherent. Import population data, fit, and
    extract the gate fidelity.
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords
        self.PP_2D = np.array([])  # purity probability
        self.NoF = 1    # number of files

        self.params_2D = None
        self.cov_2D = None
        self.std_2D = None

        self.r_inc = None
        self.r_inc_std = None

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             # data_type1='Purity_Probability'
                             data_type1='Purity_Normalized_Probability'
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        for i in range(len(self.NoC)):
            self.NoC[i] = self.NoC[i] - 1
        PP_2D = data0[data_type1]   # purity probability

        for f in file_path[1:]:
            self.NoF = self.NoF + 1
            data = noiselib.loadmat(f)
            PP_1D = data[data_type1]
            PP_2D = np.vstack((PP_2D, PP_1D))

        self.PP_2D = PP_2D

    def data_analysis(self):
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        PP_2D = self.PP_2D
        NoC = self.NoC

        params_2D = []
        cov_2D = []
        std_2D = []

        p0 = [0.9, 0.985, 0.01]
        if self.NoF == 1:
            PP_1D = PP_2D
            params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
            params_2D.append(params)
            cov_2D.append(cov)
            std_2D.append(np.sqrt(np.diag(cov)))
        else:
            for i in range(len(PP_2D)):
                PP_1D = PP_2D[i]
                params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
                params_2D.append(params)
                cov_2D.append(cov)
                std_2D.append(np.sqrt(np.diag(cov)))
        ### 2D data update
        params_2D = np.asarray(params_2D)
        cov_2D = np.asarray(cov_2D)
        self.params_2D = params_2D
        self.cov_2D = cov_2D
        self.std_2D = std_2D

    def _function(self, NoC, A, u, B):
        return A * u ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the incoherent error and std or se
        :return:
        """
        params_2D = self.params_2D
        std_2D = self.std_2D
        print('std_2D=', std_2D)
        if self.NoF == 1:
            u = params_2D[0][1]
            r_inc = 0.5*(1-np.sqrt(u))
            r_inc_std = r_inc * 0.5 * std_2D[0][1]/u
        else:
            u_list = params_2D[:, 1]
            r = [0.5*(1-np.sqrt(ur)) for ur in u_list]
            print('r=', r)
            r_inc = np.mean(r)
            r_inc_std = np.std(r)/np.sqrt(len(r))

        self.r_inc = r_inc
        self.r_inc_std = r_inc_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)
        PP_2D = self.PP_2D
        params = self.params_2D
        std_2D = self.std_2D
        # print('params=', params)
        # print('std_2D=', std_2D)

        r_inc = self.r_inc
        r_inc_std = self.r_inc_std

        PP_2D_mean = np.mean(PP_2D, axis=0)
        PP_2D_se = np.std(PP_2D, axis=0)/len(PP_2D)
        params_mean = np.mean(params, axis=0)

        mpl.rc('font', family='Arial')
        plt.figure(figsize=(8, 5))
        label_font = 20
        tick_font = 20
        legend_font = 16


        plt.errorbar(NoC, PP_2D_mean, yerr=PP_2D_se,
                     fmt="o", color='k', capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *params_mean), 'k')
        print('r_inc=', r_inc)
        print('r_inc_std=', r_inc_std)

        """real plot"""

        # plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Purity Probability')
        plt.xlim([0, 100])
        plt.ylim([0.0, 0.8])
        plt.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
        plt.xlabel('m - Number of Cliffords', fontsize=label_font)
        plt.ylabel('Sequence Purity', fontsize=label_font)
        path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
        # plt.savefig(path + '\Purity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()


class Purity_Paper_ErrorBudget(object):
    """
    Analyze the qubit gate incoherent. Import population data, fit, and
    extract the gate fidelity.
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords
        self.PP_2D = np.array([])  # purity probability
        self.NoF = 1    # number of files

        self.params_2D = None
        self.cov_2D = None
        self.std_2D = None

        self.r_inc = None
        self.r_inc_std = None

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             # data_type1='Purity_Probability'
                             data_type1='Purity_Normalized_Probability'
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        for i in range(len(self.NoC)):
            self.NoC[i] = self.NoC[i] - 1
        PP_2D = data0[data_type1]   # purity probability

        for f in file_path[1:]:
            self.NoF = self.NoF + 1
            data = noiselib.loadmat(f)
            PP_1D = data[data_type1]
            PP_2D = np.vstack((PP_2D, PP_1D))

        self.PP_2D = PP_2D

    def data_analysis(self):
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        PP_2D = self.PP_2D
        NoC = self.NoC

        params_2D = []
        cov_2D = []
        std_2D = []

        p0 = [0.9, 0.985, 0.01]
        if self.NoF == 1:
            PP_1D = PP_2D
            params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
            params_2D.append(params)
            cov_2D.append(cov)
            std_2D.append(np.sqrt(np.diag(cov)))
        else:
            for i in range(len(PP_2D)):
                PP_1D = PP_2D[i]
                params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
                params_2D.append(params)
                cov_2D.append(cov)
                std_2D.append(np.sqrt(np.diag(cov)))
        ### 2D data update
        params_2D = np.asarray(params_2D)
        cov_2D = np.asarray(cov_2D)
        self.params_2D = params_2D
        self.cov_2D = cov_2D
        self.std_2D = std_2D

    def _function(self, NoC, A, u, B):
        return A * u ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the incoherent error and std or se
        :return:
        """
        params_2D = self.params_2D
        std_2D = self.std_2D
        print('std_2D=', std_2D)
        if self.NoF == 1:
            u = params_2D[0][1]
            r_inc = 0.5*(1-np.sqrt(u))
            r_inc_std = r_inc * 0.5 * std_2D[0][1]/u
        else:
            u_list = params_2D[:, 1]
            r = [0.5*(1-np.sqrt(ur)) for ur in u_list]
            print('r=', r)
            r_inc = np.mean(r)
            r_inc_std = np.std(r)/np.sqrt(len(r))

        self.r_inc = r_inc
        self.r_inc_std = r_inc_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)
        PP_2D = self.PP_2D
        params = self.params_2D
        std_2D = self.std_2D
        # print('params=', params)
        # print('std_2D=', std_2D)

        r_inc = self.r_inc
        r_inc_std = self.r_inc_std

        PP_2D_mean = np.mean(PP_2D, axis=0)
        PP_2D_se = np.std(PP_2D, axis=0)/len(PP_2D)
        params_mean = np.mean(params, axis=0)

        print('r_inc=', r_inc)
        print('r_inc_std=', r_inc_std)

        mpl.rc('font', family='Arial')
        # fig, axs = plt.subplots(2, figsize=(8, 8),
        #                         gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.15})
        fig, axs = plt.subplots(2, figsize=(8, 10))
        label_font = 20
        tick_font = 20
        legend_font = 14


        axs[0].errorbar(NoC, PP_2D_mean, yerr=PP_2D_se,
                     fmt="o", color='k', capsize=3.0, ms=5.0)
        axs[0].plot(NoC_fit, self._function(NoC_fit, *params_mean), 'k')


        axs[0].set_xlim([0, 100])
        axs[0].set_ylim([0.0, 0.8])
        axs[0].tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
        axs[0].set_xlabel('m - Number of Cliffords', fontsize=label_font)
        axs[0].set_ylabel('Sequence Purity', fontsize=label_font)

        time = np.arange(0, 101, 1)
        t_Cliff = 90.4
        T1_base = 26.0  # us
        # Twhite_base = 2 * T1_base  # us
        Twhite_base = 20.0  # us
        Texp_base = 3 / (1 / T1_base + 1 / Twhite_base)
        k_base = 1 / (Texp_base * 1000) * 100  # from Zijun Chen's thesis

        r_incoherent = 0.964  # 0.964% per clifford gate from purity

        # axs[1].plot(time, time * r_incoherent / t_Cliff, 'k-',
        #          label='Purity incoherent error= {:.3f} % (per 10 ns)'
        #          .format(r_incoherent / t_Cliff * 10))

        # axs[1].plot(time, time * k_base, 'k--',
        #          label='Base incoherent error= {:.3f} % (per 10 ns)'
        #          .format(k_base * 10))
        plt.errorbar(t_Cliff, 100-98.810, yerr=0.09, fmt="o", color='r', capsize=3.0, ms=8, label='Total Clifford error')
        plt.errorbar(t_Cliff, 0.96429, yerr=0.0156259, fmt="d", color='k', capsize=3.0, ms=10, label='Incoherent error')
        # plt.errorbar(t_Cliff, t_Cliff * k_base, fmt="*", color='b', capsize=3.0, ms=12, label='Cliff base incoherent')
        plt.errorbar(time, time * k_base, fmt="--", color='b', capsize=3.0, ms=12, label='Baseline incoherent error')

        plt.errorbar(80, 100-99.151, yerr=0.169, fmt="s", color='r', capsize=3.0, ms=8, label='X')
        plt.errorbar(40, 100-99.503, yerr=0.160, fmt="s", color='y', capsize=3.0, ms=8, label='X/2')
        plt.errorbar(40, 100-99.383, yerr=0.128, fmt="s", color='g', capsize=3.0, ms=8, label='-X/2')

        plt.errorbar(80, 100-99.125, yerr=0.179, fmt="v", color='b', capsize=3.0, ms=8, label='Y')
        plt.errorbar(40, 100-99.377, yerr=0.180, fmt="v", color='c', capsize=3.0, ms=8, label='Y/2')
        plt.errorbar(40, 100-99.338, yerr=0.169, fmt="v", color='m', capsize=3.0, ms=8, label='-Y/2')


        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        # remove the errorbars
        handles = [h[0] for h in handles]
        plt.legend(handles, labels, loc=2, prop={'size': legend_font}, frameon=False)
        axs[1].set_xlim([0, 100])
        axs[1].set_ylim([0.0, 1.3])
        axs[1].set_yticks([0.0, 0.4, 0.8, 1.2])
        axs[1].set_xlabel('Gate Time (ns)', fontsize=label_font)
        axs[1].set_ylabel(r'Gate Error ($\times 10^{-2}$)', fontsize=label_font)
        axs[1].tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)


        path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
        plt.savefig(path + '\Purity.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        plt.show()


class Purity(object):
    """
    Analyze the qubit gate incoherent. Import population data, fit, and
    extract the gate fidelity.
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords
        self.PP_2D = np.array([])  # purity probability
        self.NoF = 1    # number of files

        self.params_2D = None
        self.cov_2D = None
        self.std_2D = None

        self.r_inc = None
        self.r_inc_std = None

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             # data_type1='Purity_Probability'
                             data_type1='Purity_Normalized_Probability'
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        for i in range(len(self.NoC)):
            self.NoC[i] = self.NoC[i] - 1
        # PP_2D = data0[data_type1]   # purity probability
        PP_2D = data0[data_type1]   # purity probability
        # print('PP_2D=', PP_2D)
        # for i in range(len(PP_2D)):
        #     PP_2D[i] = ((PP_2D[i]+1.0)**2.0)/4.0
        # PP_2D = np.array([(y+1)**2.0/4] for y in PP_2D)
        # print('PP_2D=', PP_2D)

        for f in file_path[1:]:
            self.NoF = self.NoF + 1
            data = noiselib.loadmat(f)
            # PP_1D = np.array(data[data_type1])
            PP_1D = data[data_type1]
            # for i in range(len(PP_1D)):
            #     PP_1D[i] = ((PP_1D[i] + 1.0) ** 2.0) / 4.0
            PP_2D = np.vstack((PP_2D, PP_1D))

        self.PP_2D = PP_2D

    def data_analysis(self):
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        PP_2D = self.PP_2D
        NoC = self.NoC

        params_2D = []
        cov_2D = []
        std_2D = []

        p0 = [0.9, 0.985, 0.01]
        if self.NoF == 1:
            PP_1D = PP_2D
            params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
            params_2D.append(params)
            cov_2D.append(cov)
            std_2D.append(np.sqrt(np.diag(cov)))
        else:
            for i in range(len(PP_2D)):
                PP_1D = PP_2D[i]
                params, cov = curve_fit(self._function, NoC, PP_1D, p0=p0)
                params_2D.append(params)
                cov_2D.append(cov)
                std_2D.append(np.sqrt(np.diag(cov)))
        ### 2D data update
        params_2D = np.asarray(params_2D)
        cov_2D = np.asarray(cov_2D)
        self.params_2D = params_2D
        self.cov_2D = cov_2D
        self.std_2D = std_2D

    def _function(self, NoC, A, u, B):
        # print('NoC=', NoC)
        # for i in range(len(NoC)):
        #     NoC[i] = NoC[i] - 1
        # print('NoC=', NoC)
        return A * u ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the incoherent error and std or se
        :return:
        """
        params_2D = self.params_2D
        std_2D = self.std_2D
        print('std_2D=', std_2D)
        if self.NoF == 1:
            u = params_2D[0][1]
            r_inc = 0.5*(1-np.sqrt(u))
            r_inc_std = r_inc * 0.5 * std_2D[0][1]/u
        else:
            u_list = params_2D[:, 1]
            r = [0.5*(1-np.sqrt(ur)) for ur in u_list]
            print('r=', r)
            r_inc = np.mean(r)
            r_inc_std = np.std(r)/np.sqrt(len(r))

        self.r_inc = r_inc
        self.r_inc_std = r_inc_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)
        PP_2D = self.PP_2D
        params = self.params_2D
        std_2D = self.std_2D
        print('params=', params)
        print('std_2D=', std_2D)

        r_inc = self.r_inc
        r_inc_std = self.r_inc_std

        """real plot"""
        if self.NoF == 1:
            plt.errorbar(NoC, PP_2D, fmt="o")
            plt.plot(NoC_fit, self._function(NoC_fit, *params[0]), 'k')
        else:
            for i in range(self.NoF):
                plt.errorbar(NoC, PP_2D[i])
                plt.plot(NoC_fit, self._function(NoC_fit, *params[i]))
        # plt.figure()
        # plt.errorbar(
        #     NoC, NIP_avg, yerr=NIP_ste, fmt="o", color='k', capsize=3.0,
        #     ms=5.0,
        #     label='Interleaved Gate: {}, AvgError = {:.3f} %'.format('None',
        #                                                              r_NI * 100))
        # plt.plot(NoC_fit, self._function(NoC_fit, *params), 'k')
        #
        # plt.errorbar(
        #     NoC, IP_avg, yerr=IP_ste, fmt="o", color='r', capsize=3.0, ms=5.0,
        #     label='Interleaved Gate: {}, AvgError = {:.3f} %'.format(IntGate,
        #                                                              r_I * 100))
        # plt.plot(NoC_fit, self._function(NoC_fit, *params_IP), 'r')

        plt.title('Incoherent error = {:.3f} $\pm$ {:.3f} %'.
                  format(r_inc * 100, r_inc_std * 100))

        plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Purity Probability')
        plt.show()


class RB_AllGates_Paper(object):
    """
    Analyze the qubit gate fidelity. Import population data, fit, and
    extract the gate fidelity. This is for all gates
    ['X', 'X/2', '-X/2', 'Y', 'Y/2', '-Y/2'] with one ref sequence
    For paper plot
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords

        self.Ref_sf_2D = np.array([])  # 2D ref sequence fidelity data
        self.Ref_sf_avg = np.array([])  # ref sequence fidelity avg
        self.Ref_sf_std = np.array(
            [])  # ref sequence fidelity standard deviation
        self.Ref_sf_ste = np.array([])  # ref sequence fidelity standard error

        self.X_sf_2D = np.array([])
        self.Y_sf_2D = np.array([])
        self.XOver2_sf_2D = np.array([])
        self.XOver2Neg_sf_2D = np.array([])
        self.YOver2_sf_2D = np.array([])
        self.YOver2Neg_sf_2D = np.array([])

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             data_type1='Ref_Sequence_Fidelity',
                             data_type2='X_Sequence_Fidelity',
                             data_type3='Y_Sequence_Fidelity',
                             data_type4='XOver2_Sequence_Fidelity',
                             data_type5='XOver2Neg_Sequence_Fidelity',
                             data_type6='YOver2_Sequence_Fidelity',
                             data_type7='YOver2Neg_Sequence_Fidelity',
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        Ref_sf_2D = data0[data_type1]
        X_sf_2D = data0[data_type2]
        Y_sf_2D = data0[data_type3]
        XOver2_sf_2D = data0[data_type3]
        XOver2Neg_sf_2D = data0[data_type3]
        YOver2_sf_2D = data0[data_type3]
        YOver2Neg_sf_2D = data0[data_type3]

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            Ref_sf_1D = np.array(data[data_type1])
            Ref_sf_2D = np.vstack((Ref_sf_2D, Ref_sf_1D))
            X_sf_1D = np.array(data[data_type2])
            X_sf_2D = np.vstack((X_sf_2D, X_sf_1D))
            Y_sf_1D = np.array(data[data_type3])
            Y_sf_2D = np.vstack((Y_sf_2D, Y_sf_1D))
            XOver2_sf_1D = np.array(data[data_type4])
            XOver2_sf_2D = np.vstack((XOver2_sf_2D, XOver2_sf_1D))
            XOver2Neg_sf_1D = np.array(data[data_type5])
            XOver2Neg_sf_2D = np.vstack((XOver2Neg_sf_2D, XOver2Neg_sf_1D))
            YOver2_sf_1D = np.array(data[data_type6])
            YOver2_sf_2D = np.vstack((YOver2_sf_2D, YOver2_sf_1D))
            YOver2Neg_sf_1D = np.array(data[data_type7])
            YOver2Neg_sf_2D = np.vstack((YOver2Neg_sf_2D, YOver2Neg_sf_1D))

        self.Ref_sf_2D = Ref_sf_2D
        self.X_sf_2D = X_sf_2D
        self.Y_sf_2D = Y_sf_2D
        self.XOver2_sf_2D = XOver2_sf_2D
        self.XOver2Neg_sf_2D = XOver2Neg_sf_2D
        self.YOver2_sf_2D = YOver2_sf_2D
        self.YOver2Neg_sf_2D = YOver2Neg_sf_2D

    def data_analysis(self):
        self._convert_data()
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _convert_data(self):
        """
        Make 2D raw data to 1D average
        :return:
        """
        self.NoC = self.NoC
        m = len(self.NoC)
        NoS = len(self.Ref_sf_2D)  # number of seqs

        self.Ref_sf_avg = np.average(self.Ref_sf_2D, axis=0)
        self.Ref_sf_std = np.std(self.Ref_sf_2D, axis=0)
        self.Ref_sf_ste = self.Ref_sf_std / np.sqrt(NoS)

        self.X_sf_avg = np.average(self.X_sf_2D, axis=0)
        self.X_sf_std = np.std(self.X_sf_2D, axis=0)
        self.X_sf_ste = self.X_sf_std / np.sqrt(NoS)

        self.Y_sf_avg = np.average(self.Y_sf_2D, axis=0)
        self.Y_sf_std = np.std(self.Y_sf_2D, axis=0)
        self.Y_sf_ste = self.Y_sf_std / np.sqrt(NoS)

        self.XOver2_sf_avg = np.average(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_std = np.std(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_ste = self.XOver2_sf_std / np.sqrt(NoS)

        self.XOver2Neg_sf_avg = np.average(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_std = np.std(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_ste = self.XOver2Neg_sf_std / np.sqrt(NoS)

        self.YOver2_sf_avg = np.average(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_std = np.std(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_ste = self.YOver2_sf_std / np.sqrt(NoS)

        self.YOver2Neg_sf_avg = np.average(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_std = np.std(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_ste = self.YOver2Neg_sf_std / np.sqrt(NoS)

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        NoC = self.NoC
        p0 = [0.5, 0.985, 0.5]

        Ref_sf_avg = self.Ref_sf_avg
        params_Ref, covariance_Ref = curve_fit(self._function, NoC, Ref_sf_avg,
                                               p0=p0)
        self.params_Ref = params_Ref
        self.err_Ref = np.sqrt(np.diag(covariance_Ref))

        X_sf_avg = self.X_sf_avg
        params_X, covariance_X = curve_fit(self._function, NoC, X_sf_avg,
                                           p0=p0)
        self.params_X = params_X
        self.err_X = np.sqrt(np.diag(covariance_X))

        Y_sf_avg = self.Y_sf_avg
        params_Y, covariance_Y = curve_fit(self._function, NoC, Y_sf_avg,
                                           p0=p0)
        self.params_Y = params_Y
        self.err_Y = np.sqrt(np.diag(covariance_Y))

        XOver2_sf_avg = self.XOver2_sf_avg
        params_XOver2, covariance_XOver2 = curve_fit(self._function, NoC,
                                                     XOver2_sf_avg, p0=p0)
        self.params_XOver2 = params_XOver2
        self.err_XOver2 = np.sqrt(np.diag(covariance_XOver2))

        XOver2Neg_sf_avg = self.XOver2Neg_sf_avg
        params_XOver2Neg, covariance_XOver2Neg = curve_fit(self._function, NoC,
                                                           XOver2Neg_sf_avg,
                                                           p0=p0)
        self.params_XOver2Neg = params_XOver2Neg
        self.err_XOver2Neg = np.sqrt(np.diag(covariance_XOver2Neg))

        YOver2_sf_avg = self.YOver2_sf_avg
        params_YOver2, covariance_YOver2 = curve_fit(self._function, NoC,
                                                     YOver2_sf_avg, p0=p0)
        self.params_YOver2 = params_YOver2
        self.err_YOver2 = np.sqrt(np.diag(covariance_YOver2))

        YOver2Neg_sf_avg = self.YOver2Neg_sf_avg
        params_YOver2Neg, covariance_YOver2Neg = curve_fit(self._function, NoC,
                                                           YOver2Neg_sf_avg,
                                                           p0=p0)
        self.params_YOver2Neg = params_YOver2Neg
        self.err_YOver2Neg = np.sqrt(np.diag(covariance_YOver2Neg))

    def _function(self, NoC, A, p, B):
        return A * p ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the (non) interleaved seq error and fidelity, and
        find out the interleaved gate's fidelity
        :return:
        """
        p_Ref = self.params_Ref[1]
        err_Ref = self.err_Ref[1]
        r_Ref = (1 - p_Ref) / 2  # non-interleaved sequence error
        F_Ref = 1 - r_Ref
        F_Ref_std = err_Ref
        self.r_Ref = r_Ref  # avg gate error
        self.F_Ref = F_Ref  # avg gate fidelity
        self.F_Ref_std = F_Ref_std  # avg gate fidelity standard error

        p_X = self.params_X[1]
        err_X = self.err_X[1]
        r_X = (1 - p_X) / 2  # interleaved sequence error
        F_X = (1 + p_X / p_Ref) / 2.0  # X gate fidelity
        F_X_std = F_X * np.sqrt((err_Ref / p_Ref) ** 2 + (err_X / p_X) ** 2)
        self.r_X = r_X  # avg gate+X error
        self.F_X = F_X  # X gate fidelity
        self.F_X_std = F_X_std  # X gate fidelity std

        p_Y = self.params_Y[1]
        err_Y = self.err_Y[1]
        r_Y = (1 - p_Y) / 2  # interleaved sequence error
        F_Y = (1 + p_Y / p_Ref) / 2.0  # X gate fidelity
        F_Y_std = F_Y * np.sqrt((err_Ref / p_Ref) ** 2 + (err_Y / p_Y) ** 2)
        self.r_Y = r_Y  # avg gate+X error
        self.F_Y = F_Y  # X gate fidelity
        self.F_Y_std = F_Y_std  # X gate fidelity std

        p_XOver2 = self.params_XOver2[1]
        err_XOver2 = self.err_XOver2[1]
        r_XOver2 = (1 - p_XOver2) / 2
        F_XOver2 = (1 + p_XOver2 / p_Ref) / 2.0
        F_XOver2_std = F_XOver2 * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_XOver2 / p_XOver2) ** 2)
        self.r_XOver2 = r_XOver2
        self.F_XOver2 = F_XOver2
        self.F_XOver2_std = F_XOver2_std

        p_XOver2Neg = self.params_XOver2Neg[1]
        err_XOver2Neg = self.err_XOver2Neg[1]
        r_XOver2Neg = (1 - p_XOver2Neg) / 2
        F_XOver2Neg = (1 + p_XOver2Neg / p_Ref) / 2.0
        F_XOver2Neg_std = F_XOver2Neg * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_XOver2Neg / p_XOver2Neg) ** 2)
        self.r_XOver2Neg = r_XOver2Neg
        self.F_XOver2Neg = F_XOver2Neg
        self.F_XOver2Neg_std = F_XOver2Neg_std

        p_YOver2 = self.params_YOver2[1]
        err_YOver2 = self.err_YOver2[1]
        r_YOver2 = (1 - p_YOver2) / 2
        F_YOver2 = (1 + p_YOver2 / p_Ref) / 2.0
        F_YOver2_std = F_YOver2 * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2 / p_YOver2) ** 2)
        self.r_YOver2 = r_YOver2
        self.F_YOver2 = F_YOver2
        self.F_YOver2_std = F_YOver2_std

        p_YOver2Neg = self.params_YOver2Neg[1]
        err_YOver2Neg = self.err_YOver2Neg[1]
        r_YOver2Neg = (1 - p_YOver2Neg) / 2
        F_YOver2Neg = (1 + p_YOver2Neg / p_Ref) / 2.0
        F_YOver2Neg_std = F_YOver2Neg * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2Neg / p_YOver2Neg) ** 2)
        self.r_YOver2Neg = r_YOver2Neg
        self.F_YOver2Neg = F_YOver2Neg
        self.F_YOver2Neg_std = F_YOver2Neg_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(2, 97, 1)

        """real plot"""

        # rcParams['mathtext.fontset'] = 'custom'
        # rcParams['mathtext.it'] = 'Arial:italic'
        # rcParams['mathtext.rm'] = 'Arial'
        # mat.rc('font', family='Times New Roman')
        mpl.rc('font', family='Arial')

        plt.figure(figsize=(8, 5))

        # plt.errorbar(0, 0, linestyle='None', ms=0, label='Gate')

        plt.errorbar(NoC, self.Ref_sf_avg, yerr=self.Ref_sf_ste, fmt="o",
                     color='k', capsize=3.0, ms=5.0, label='Ref')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Ref), 'k')

        plt.errorbar(NoC, self.X_sf_avg, yerr=self.X_sf_ste, fmt="s", color='r'
                     , capsize=3.0, ms=5.0, label='X')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_X), 'r')

        plt.errorbar(NoC, self.XOver2_sf_avg, yerr=self.XOver2_sf_ste, fmt="s",
                     color='y', capsize=3.0, ms=5.0, label='X/2')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2), 'y')

        plt.errorbar(NoC, self.XOver2Neg_sf_avg, yerr=self.XOver2Neg_sf_ste,
                     fmt="s", color='g', capsize=3.0, ms=5.0,  label='-X/2')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2Neg), 'g')

        plt.errorbar(NoC, self.Y_sf_avg, yerr=self.Y_sf_ste, fmt="v", color='b'
                     , capsize=3.0, ms=5.0, label='Y')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Y), 'b')

        plt.errorbar(NoC, self.YOver2_sf_avg, yerr=self.YOver2_sf_ste, fmt="v",
                     color='c', capsize=3.0, ms=5.0, label='Y/2')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2), 'c')

        plt.errorbar(NoC, self.YOver2Neg_sf_avg, yerr=self.YOver2Neg_sf_ste,
                     fmt="v", color='m', capsize=3.0, ms=5.0, label='-Y/2')
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2Neg), 'm')

        # plt.errorbar(0, 0, linestyle='None', ms=0, label='Fidelity')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.988(1)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.992(2)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.995(2)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.994(1)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.991(2)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.994(2)')
        # plt.errorbar(0, 0, linestyle='None', ms=0, label='0.993(2)')

        label_font = 20
        tick_font = 20
        legend_font = 16

        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        # remove the errorbars
        handles = [h[0] for h in handles]

        ax.legend(handles, labels, loc=(0.64, 0.32), frameon=False, prop={'size': legend_font}, handletextpad=0)
        plt.xlim([0, 100])
        # plt.ylim([0.33, 0.9]) # over4
        plt.ylim([0.31, 0.94])   # over2
        plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        plt.tick_params(labelsize=tick_font, direction="in", width=1.5, length=6)
        plt.xlabel('m - Number of Cliffords', fontsize=label_font)
        plt.ylabel('Sequence Fidelity', fontsize=label_font)

        path = 'Z:\mcdermott-group\data\sfq\SFQMCMPaperWriting\FromPython'
        # plt.savefig(path + '\RBIRB.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        # plt.savefig(path + '\RBIRBOver2.pdf', format='pdf', bbox_inches='tight', dpi=1200)
        # plt.show()


class RB_AllGates(object):
    """
    Analyze the qubit gate fidelity. Import population data, fit, and
    extract the gate fidelity. This is for all gates
    ['X', 'X/2', '-X/2', 'Y', 'Y/2', '-Y/2'] with one ref sequence
    """

    def __init__(self):
        self.NoC = np.array([])  # number of cliffords

        self.Ref_sf_2D = np.array([])  # 2D ref sequence fidelity data
        self.Ref_sf_avg = np.array([])  # ref sequence fidelity avg
        self.Ref_sf_std = np.array(
            [])  # ref sequence fidelity standard deviation
        self.Ref_sf_ste = np.array([])  # ref sequence fidelity standard error

        self.X_sf_2D = np.array([])
        self.Y_sf_2D = np.array([])
        self.XOver2_sf_2D = np.array([])
        self.XOver2Neg_sf_2D = np.array([])
        self.YOver2_sf_2D = np.array([])
        self.YOver2Neg_sf_2D = np.array([])

    def add_data_from_matlab(self, file_path,
                             data_type0='Number_of_Cliffords',
                             data_type1='Ref_Sequence_Fidelity',
                             data_type2='X_Sequence_Fidelity',
                             data_type3='Y_Sequence_Fidelity',
                             data_type4='XOver2_Sequence_Fidelity',
                             data_type5='XOver2Neg_Sequence_Fidelity',
                             data_type6='YOver2_Sequence_Fidelity',
                             data_type7='YOver2Neg_Sequence_Fidelity',
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.NoC = data0[data_type0]
        Ref_sf_2D = data0[data_type1]
        X_sf_2D = data0[data_type2]
        Y_sf_2D = data0[data_type3]
        XOver2_sf_2D = data0[data_type3]
        XOver2Neg_sf_2D = data0[data_type3]
        YOver2_sf_2D = data0[data_type3]
        YOver2Neg_sf_2D = data0[data_type3]

        for f in file_path[1:]:
            data = noiselib.loadmat(f)
            Ref_sf_1D = np.array(data[data_type1])
            Ref_sf_2D = np.vstack((Ref_sf_2D, Ref_sf_1D))
            X_sf_1D = np.array(data[data_type2])
            X_sf_2D = np.vstack((X_sf_2D, X_sf_1D))
            Y_sf_1D = np.array(data[data_type3])
            Y_sf_2D = np.vstack((Y_sf_2D, Y_sf_1D))
            XOver2_sf_1D = np.array(data[data_type4])
            XOver2_sf_2D = np.vstack((XOver2_sf_2D, XOver2_sf_1D))
            XOver2Neg_sf_1D = np.array(data[data_type5])
            XOver2Neg_sf_2D = np.vstack((XOver2Neg_sf_2D, XOver2Neg_sf_1D))
            YOver2_sf_1D = np.array(data[data_type6])
            YOver2_sf_2D = np.vstack((YOver2_sf_2D, YOver2_sf_1D))
            YOver2Neg_sf_1D = np.array(data[data_type7])
            YOver2Neg_sf_2D = np.vstack((YOver2Neg_sf_2D, YOver2Neg_sf_1D))

        self.Ref_sf_2D = Ref_sf_2D
        self.X_sf_2D = X_sf_2D
        self.Y_sf_2D = Y_sf_2D
        self.XOver2_sf_2D = XOver2_sf_2D
        self.XOver2Neg_sf_2D = XOver2Neg_sf_2D
        self.YOver2_sf_2D = YOver2_sf_2D
        self.YOver2Neg_sf_2D = YOver2Neg_sf_2D

    def data_analysis(self):
        self._convert_data()
        self._fit_data()
        self._extractErrorAndFidelity()
        return

    def _convert_data(self):
        """
        Make 2D raw data to 1D average
        :return:
        """
        self.NoC = self.NoC
        m = len(self.NoC)
        NoS = len(self.Ref_sf_2D)  # number of seqs

        self.Ref_sf_avg = np.average(self.Ref_sf_2D, axis=0)
        self.Ref_sf_std = np.std(self.Ref_sf_2D, axis=0)
        self.Ref_sf_ste = self.Ref_sf_std / np.sqrt(NoS)

        self.X_sf_avg = np.average(self.X_sf_2D, axis=0)
        self.X_sf_std = np.std(self.X_sf_2D, axis=0)
        self.X_sf_ste = self.X_sf_std / np.sqrt(NoS)

        self.Y_sf_avg = np.average(self.Y_sf_2D, axis=0)
        self.Y_sf_std = np.std(self.Y_sf_2D, axis=0)
        self.Y_sf_ste = self.Y_sf_std / np.sqrt(NoS)

        self.XOver2_sf_avg = np.average(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_std = np.std(self.XOver2_sf_2D, axis=0)
        self.XOver2_sf_ste = self.XOver2_sf_std / np.sqrt(NoS)

        self.XOver2Neg_sf_avg = np.average(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_std = np.std(self.XOver2Neg_sf_2D, axis=0)
        self.XOver2Neg_sf_ste = self.XOver2Neg_sf_std / np.sqrt(NoS)

        self.YOver2_sf_avg = np.average(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_std = np.std(self.YOver2_sf_2D, axis=0)
        self.YOver2_sf_ste = self.YOver2_sf_std / np.sqrt(NoS)

        self.YOver2Neg_sf_avg = np.average(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_std = np.std(self.YOver2Neg_sf_2D, axis=0)
        self.YOver2Neg_sf_ste = self.YOver2Neg_sf_std / np.sqrt(NoS)

    def _fit_data(self):
        """
        Fit the (non) interleaved sequences
        :return:
        """
        NoC = self.NoC
        p0 = [0.5, 0.985, 0.5]

        Ref_sf_avg = self.Ref_sf_avg
        params_Ref, covariance_Ref = curve_fit(self._function, NoC, Ref_sf_avg,
                                               p0=p0)
        self.params_Ref = params_Ref
        self.err_Ref = np.sqrt(np.diag(covariance_Ref))

        X_sf_avg = self.X_sf_avg
        params_X, covariance_X = curve_fit(self._function, NoC, X_sf_avg,
                                           p0=p0)
        self.params_X = params_X
        self.err_X = np.sqrt(np.diag(covariance_X))

        Y_sf_avg = self.Y_sf_avg
        params_Y, covariance_Y = curve_fit(self._function, NoC, Y_sf_avg,
                                           p0=p0)
        self.params_Y = params_Y
        self.err_Y = np.sqrt(np.diag(covariance_Y))

        XOver2_sf_avg = self.XOver2_sf_avg
        params_XOver2, covariance_XOver2 = curve_fit(self._function, NoC,
                                                     XOver2_sf_avg, p0=p0)
        self.params_XOver2 = params_XOver2
        self.err_XOver2 = np.sqrt(np.diag(covariance_XOver2))

        XOver2Neg_sf_avg = self.XOver2Neg_sf_avg
        params_XOver2Neg, covariance_XOver2Neg = curve_fit(self._function, NoC,
                                                           XOver2Neg_sf_avg,
                                                           p0=p0)
        self.params_XOver2Neg = params_XOver2Neg
        self.err_XOver2Neg = np.sqrt(np.diag(covariance_XOver2Neg))

        YOver2_sf_avg = self.YOver2_sf_avg
        params_YOver2, covariance_YOver2 = curve_fit(self._function, NoC,
                                                     YOver2_sf_avg, p0=p0)
        self.params_YOver2 = params_YOver2
        self.err_YOver2 = np.sqrt(np.diag(covariance_YOver2))

        YOver2Neg_sf_avg = self.YOver2Neg_sf_avg
        params_YOver2Neg, covariance_YOver2Neg = curve_fit(self._function, NoC,
                                                           YOver2Neg_sf_avg,
                                                           p0=p0)
        self.params_YOver2Neg = params_YOver2Neg
        self.err_YOver2Neg = np.sqrt(np.diag(covariance_YOver2Neg))

    def _function(self, NoC, A, p, B):
        return A * p ** NoC + B

    def _extractErrorAndFidelity(self):
        """
        To extract the (non) interleaved seq error and fidelity, and
        find out the interleaved gate's fidelity
        :return:
        """
        p_Ref = self.params_Ref[1]
        err_Ref = self.err_Ref[1]
        r_Ref = (1 - p_Ref) / 2  # non-interleaved sequence error
        F_Ref = 1 - r_Ref
        F_Ref_std = err_Ref
        self.r_Ref = r_Ref  # avg gate error
        self.F_Ref = F_Ref  # avg gate fidelity
        self.F_Ref_std = F_Ref_std  # avg gate fidelity standard error

        p_X = self.params_X[1]
        err_X = self.err_X[1]
        r_X = (1 - p_X) / 2  # interleaved sequence error
        F_X = (1 + p_X / p_Ref) / 2.0  # X gate fidelity
        F_X_std = F_X * np.sqrt((err_Ref / p_Ref) ** 2 + (err_X / p_X) ** 2)
        self.r_X = r_X  # avg gate+X error
        self.F_X = F_X  # X gate fidelity
        self.F_X_std = F_X_std  # X gate fidelity std

        p_Y = self.params_Y[1]
        err_Y = self.err_Y[1]
        r_Y = (1 - p_Y) / 2  # interleaved sequence error
        F_Y = (1 + p_Y / p_Ref) / 2.0  # X gate fidelity
        F_Y_std = F_Y * np.sqrt((err_Ref / p_Ref) ** 2 + (err_Y / p_Y) ** 2)
        self.r_Y = r_Y  # avg gate+X error
        self.F_Y = F_Y  # X gate fidelity
        self.F_Y_std = F_Y_std  # X gate fidelity std

        p_XOver2 = self.params_XOver2[1]
        err_XOver2 = self.err_XOver2[1]
        r_XOver2 = (1 - p_XOver2) / 2
        F_XOver2 = (1 + p_XOver2 / p_Ref) / 2.0
        F_XOver2_std = F_XOver2 * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_XOver2 / p_XOver2) ** 2)
        self.r_XOver2 = r_XOver2
        self.F_XOver2 = F_XOver2
        self.F_XOver2_std = F_XOver2_std

        p_XOver2Neg = self.params_XOver2Neg[1]
        err_XOver2Neg = self.err_XOver2Neg[1]
        r_XOver2Neg = (1 - p_XOver2Neg) / 2
        F_XOver2Neg = (1 + p_XOver2Neg / p_Ref) / 2.0
        F_XOver2Neg_std = F_XOver2Neg * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_XOver2Neg / p_XOver2Neg) ** 2)
        self.r_XOver2Neg = r_XOver2Neg
        self.F_XOver2Neg = F_XOver2Neg
        self.F_XOver2Neg_std = F_XOver2Neg_std

        p_YOver2 = self.params_YOver2[1]
        err_YOver2 = self.err_YOver2[1]
        r_YOver2 = (1 - p_YOver2) / 2
        F_YOver2 = (1 + p_YOver2 / p_Ref) / 2.0
        F_YOver2_std = F_YOver2 * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2 / p_YOver2) ** 2)
        self.r_YOver2 = r_YOver2
        self.F_YOver2 = F_YOver2
        self.F_YOver2_std = F_YOver2_std

        p_YOver2Neg = self.params_YOver2Neg[1]
        err_YOver2Neg = self.err_YOver2Neg[1]
        r_YOver2Neg = (1 - p_YOver2Neg) / 2
        F_YOver2Neg = (1 + p_YOver2Neg / p_Ref) / 2.0
        F_YOver2Neg_std = F_YOver2Neg * np.sqrt(
            (err_Ref / p_Ref) ** 2 + (err_YOver2Neg / p_YOver2Neg) ** 2)
        self.r_YOver2Neg = r_YOver2Neg
        self.F_YOver2Neg = F_YOver2Neg
        self.F_YOver2Neg_std = F_YOver2Neg_std

    def plot(self):
        """Data, avg and err"""
        NoC = self.NoC
        NoC_fit = np.arange(1, 97, 1)

        """real plot"""
        plt.figure()
        # plt.errorbar(
        #     NoC, NIP_avg, yerr=NIP_ste, fmt="o", color='k', capsize=3.0, ms=5.0,
        #     label='Interleaved Gate: {}, AvgError = {:.3f} %'.format('None', r_NI * 100))
        plt.errorbar(NoC, self.Ref_sf_avg, yerr=self.Ref_sf_ste, fmt="o",
                     color='k'
                     , capsize=3.0, ms=5.0)
        # plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Ref), 'k',
        #          label='Ref, AvgError = {:.3f} %'.format(self.r_Ref * 100))
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Ref), 'k',
                 label='Ref, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_Ref * 100, self.F_Ref * 100,
                         self.F_Ref_std * 100))

        plt.errorbar(NoC, self.X_sf_avg, yerr=self.X_sf_ste, fmt="o", color='r'
                     , capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_X), 'r',
                 label='X, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_X * 100, self.F_X * 100, self.F_X_std * 100))
        #
        plt.errorbar(NoC, self.Y_sf_avg, yerr=self.Y_sf_ste, fmt="o", color='b'
                     , capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_Y), 'b',
                 label='Y, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_Y * 100, self.F_Y * 100, self.F_Y_std * 100))

        plt.errorbar(NoC, self.XOver2_sf_avg, yerr=self.XOver2_sf_ste, fmt="o",
                     color='y'
                     , capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2), 'y',
                 label='X/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_XOver2 * 100, self.F_XOver2 * 100,
                         self.F_XOver2_std * 100))

        plt.errorbar(NoC, self.XOver2Neg_sf_avg, yerr=self.XOver2Neg_sf_ste,
                     fmt="o", color='g', capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_XOver2Neg), 'g',
                 label='-X/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_XOver2Neg * 100, self.F_XOver2Neg * 100,
                         self.F_XOver2Neg_std * 100))

        plt.errorbar(NoC, self.YOver2_sf_avg, yerr=self.YOver2_sf_ste, fmt="o",
                     color='b', capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2), 'b',
                 label='Y/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_YOver2 * 100, self.F_YOver2 * 100,
                         self.F_YOver2_std * 100))

        plt.errorbar(NoC, self.YOver2Neg_sf_avg, yerr=self.YOver2Neg_sf_ste,
                     fmt="o", color='r', capsize=3.0, ms=5.0)
        plt.plot(NoC_fit, self._function(NoC_fit, *self.params_YOver2Neg), 'r',
                 label='-Y/2, AvgError = {:.3f} %. Fidelity= {:.3f} $\pm$ {:.3f} %'
                 .format(self.r_YOver2Neg * 100, self.F_YOver2Neg * 100,
                         self.F_YOver2Neg_std * 100))

        plt.title('SFQ-Based Gate Fidelities')
        plt.legend()
        plt.xlabel('Number of Cliffords')
        plt.ylabel('Sequence Fidelity')
        plt.show()


class RB_ORBIT(object):
    """
    Use OBRIT to optimize parameters
    """

    def __init__(self):
        self.paramOfInterest = np.array([])  # the parameter of interest
        self.seqFidelity_2D = np.array([])  # sequence fidelity 2D
        self.seqFidelity_avg = np.array([])  # sequence fidelity avg
        self.seqFidelity_se = np.array([])  # sequence fidelity standard error
        self.NoF = 1    # number of files
        self.paramOfInterest_Name = ''

    def add_data_from_matlab(self, file_path,
                             data_type0='SFQ_Y_Gate_Phase_Offset',
                             data_type1='Sequence_Fidelity'
                             ):
        f0 = file_path[0]  # first file
        data0 = noiselib.loadmat(f0)
        self.paramOfInterest = data0[data_type0]
        self.paramOfInterest_Name = data_type0
        seqF_2D = data0[data_type1]   # sequence fidelity
        NoF = self.NoF
        for f in file_path[1:]:
            NoF = NoF + 1
            data = noiselib.loadmat(f)
            seqF_1D = data[data_type1]
            seqF_2D = np.vstack((seqF_2D, seqF_1D))

        self.seqFidelity_2D = seqF_2D
        self.NoF = NoF

        if NoF == 1:    # no avg needed
            self.seqFidelity_avg = seqF_2D
            self.seqFidelity_se = 0 * seqF_2D
        else:
            self.seqFidelity_avg = np.average(seqF_2D, axis=0)
            self.seqFidelity_se = np.std(seqF_2D, axis=0)/np.sqrt(NoF)

    def plot(self):
        paramOfInterest = self.paramOfInterest
        paramOfInterest_Name = self.paramOfInterest_Name
        seqFidelity_avg = self.seqFidelity_avg
        seqFidelity_se = self.seqFidelity_se

        plt.xlabel(paramOfInterest_Name)
        plt.ylabel('Sequence Fidelity')
        plt.errorbar(paramOfInterest, seqFidelity_avg, yerr=seqFidelity_se)
        plt.show()