import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sc
import math

class Calibration(object):
    def __init__(self, I1, Q1, I2, Q2):
        # raw IQ data for calibration
        self.I1 = I1
        self.Q1 = Q1
        self.I2 = I2
        self.Q2 = Q2
        self.states = ['g', 'e'] # state labels (used for plotting).

        # rotated IQ data
        self.angle = 0.0
        self.rI1 = None
        self.rQ1 = None
        self.rI2 = None
        self.rQ2 = None

        self.rI1mean = None
        self.rI2mean = None

        self.threshold = None

        self._get_angle()
        self._get_rot_blobs()
        self._get_threshold()

    def _get_angle(self):
        self.angle = np.arctan2(np.mean(self.Q2) - np.mean(self.Q1), np.mean(self.I1) - np.mean(self.I2))

    def _get_rot_blobs(self):
        self.rI1 = self.I1 * np.cos(self.angle) - self.Q1 * np.sin(self.angle)
        self.rQ1 = self.I1 * np.sin(self.angle) + self.Q1 * np.cos(self.angle)
        self.rI2 = self.I2 * np.cos(self.angle) - self.Q2 * np.sin(self.angle)
        self.rQ2 = self.I2 * np.sin(self.angle) + self.Q2 * np.cos(self.angle)

    def _rot_blob(self, I, Q):
        rI = I * np.cos(self.angle) - Q * np.sin(self.angle)
        rQ = I * np.sin(self.angle) + Q * np.cos(self.angle)
        return rI, rQ

    def _false_detections(self, threshold):
        '''
        From QM discriminator.py
        '''
        if np.mean(self.rI1) < np.mean(self.rI2):
            false_detections_var = np.sum(self.rI1 > threshold) + np.sum(self.rI2 < threshold)
        else:
            false_detections_var = np.sum(self.rI1 < threshold) + np.sum(self.rI2 > threshold)
        return false_detections_var

    def _get_threshold(self):
        guess = 0.5 * (np.mean(self.rI1) + np.mean(self.rI2))
        fit = sc.minimize(self._false_detections,
                          guess,
                          method='Nelder-Mead')
        self.threshold = fit.x[0]

    def _gauss(self, x, mean, sigma):
        return np.exp(-0.5 * (x - mean) ** 2 / sigma ** 2) * 1 / (sigma * np.sqrt(2 * np.pi))

    def _bimodal(self, x, mu1, sigma1, h1, mu2, sigma2, h2):
        return h1 * self._gauss(x, mu1, sigma1) + h2 * self._gauss(x, mu2, sigma2)

    def _bimodal_fix_means(self, x, sigma1, h1, sigma2, h2):
        return h1 * self._gauss(x, self.rI1mean, sigma1) + h2 * self._gauss(x, self.rI2mean, sigma2)

    def _bimodal_fix_mean1(self, x, sigma1, h1, mean2, sigma2, h2):
        return h1 * self._gauss(x, self.rI1mean, sigma1) + h2 * self._gauss(x, mean2, sigma2)

    '''
    Uppdate state labels for plotting purposes.
    '''
    def update_state_labels(self, states):
        self.states = states

    '''
    Extract population via thresholding.
    Threshold is taken to be the midpoint between the centers of the two blobs.
    '''
    def get_pop_from_threshold(self, I, Q):
        rI, rQ = self._rot_blob(I, Q)

        n1 = 0
        n2 = 0
        for i in rI:
            if np.abs(i - np.mean(self.rI1)) < np.abs(i - np.mean(self.rI2)):
                n1 += 1
            else:
                n2 += 1

        P1 = n1 / len(rI)
        P2 = n2 / len(rI)

        return P1, P2

    '''
    Extract population via thresholding. 
    Fit to threshold by minimizing false detections (a la QM two state discriminator).
    '''
    def get_pop_from_threshold_fit(self, I, Q):
        rI, rQ = self._rot_blob(I, Q)

        n1 = 0
        n2 = 0
        for i in rI:
            if np.mean(self.rI1) < self.threshold:
                if i < self.threshold:
                    n1 += 1
                else:
                    n2 += 1
            else:
                if i > self.threshold:
                    n1 += 1
                else:
                    n2 += 1

        P1 = n1 / len(rI)
        P2 = n2 / len(rI)

        return P1, P2

    '''
    Extract population by fitting to a bimodal distribution. 
    Optional: fix means.
    Closely follows Vincent's ProjectedOccupation class in calibration.py.
    '''
    def get_pop_from_fit(self, I, Q):
        rI, rQ = self._rot_blob(I, Q)
        I1_mean = np.mean(self.rI1)
        I1_std = np.std(self.rI1)
        I2_mean = np.mean(self.rI2)
        I2_std = np.std(self.rI2)

        y, x = np.histogram(rI, 200)
        x = (x[1:] + x[:-1]) / 2 # calculate centers of bins

        # estimate amplitudes
        I1_index = (np.abs(x - I1_mean)).argmin()
        I2_index = (np.abs(x - I2_mean)).argmin()
        I1_amp = np.mean(y[I1_index - 2: I1_index + 3])
        I2_amp = np.mean(y[I2_index - 2: I2_index + 3])
        guess = [I1_mean, I1_std, I1_amp, I2_mean, I2_std, I2_amp]
        print(guess)
        # put bounds on means and amplitudes
        # bounds = ([I1_mean - I1_std, -1*np.inf, 0., I2_mean - I2_std, -1*np.inf, 0.],
        #           [I1_mean + I1_std, np.inf, len(rI), I2_mean + I2_std, np.inf, len(rI)])
        # print(guess)
        # print(bounds)

        # fit
        params, cov = sc.curve_fit(self._bimodal, x, y, p0=guess,bounds=([-np.inf,0,0,-np.inf,0,0],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]),method='dogbox')
        #mu1, sigma1, h1, mu2, sigma2, h2
        P1 = (params[1] * params[2]) / (params[1] * params[2] + params[4] * params[5])
        P2 = (params[4] * params[5]) / (params[1] * params[2] + params[4] * params[5])
        # err = np.sqrt(np.diag(cov))

        return P1, P2, params, cov

    def get_pop_from_fit_fix_means(self, I, Q):
        rI, rQ = self._rot_blob(I, Q)
        I1_mean = self.rI1mean#np.mean(self.rI1)
        I1_std = np.std(self.rI1)
        I2_mean = self.rI2mean#np.mean(self.rI2)
        I2_std = np.std(self.rI2)

        y, x = np.histogram(rI, 200)
        x = (x[1:] + x[:-1]) / 2  # calculate centers of bins

        # estimate amplitudes
        I1_index = (np.abs(x - I1_mean)).argmin()
        I2_index = (np.abs(x - I2_mean)).argmin()
        I1_amp = np.mean(y[I1_index - 2: I1_index + 3])
        I2_amp = np.mean(y[I2_index - 2: I2_index + 3])

        # fit
        guess = [I1_std, I1_amp, I2_std, I2_amp]
        params, cov = sc.curve_fit(self._bimodal_fix_means, x, y, p0=guess,bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf]))
        P1 = (params[0] * params[1]) / (params[0] * params[1] + params[2] * params[3])
        P2 = (params[2] * params[3]) / (params[0] * params[1] + params[2] * params[3])

        return P1, P2, params, cov

    def get_gauss_fit(self, I, Q):
        rI, rQ = self._rot_blob(I, Q)
        y, x = np.histogram(rI, 200)
        x = (x[1:] + x[:-1]) / 2  # calculate centers of bins

        # estimate amplitudes
        rI1mean = np.mean(rI)
        I1_index = (np.abs(x - rI1mean)).argmin()
        I1_amp = np.mean(y[I1_index - 2: I1_index + 3])
        I1_std = np.std(self.rI1)

        # fit
        guess = [rI1mean,I1_std]
        params, cov = sc.curve_fit(self._gauss, x, y, p0=guess)

        return params, cov

    def get_bimodal_fit_fix_mean1(self):
        I1_mean = self.rI1mean  # np.mean(self.rI1)
        I1_std = np.std(self.rI1)
        I2_mean = np.mean(self.rI2)  # np.mean(self.rI2)
        I2_std = np.std(self.rI2)

        x1 = np.linspace(self.rI1.min(), self.rI1.max(), 200)
        x2 = np.linspace(self.rI2.min(), self.rI2.max(), 200)

        y, x = np.histogram(self.rI2, 200)
        x = (x[1:] + x[:-1]) / 2  # calculate centers of bins

        # estimate amplitudes
        I1_index = (np.abs(x - I1_mean)).argmin()
        I2_index = (np.abs(x - I2_mean)).argmin()
        I1_amp = np.mean(y[I1_index - 2: I1_index + 3])
        I2_amp = np.mean(y[I2_index - 2: I2_index + 3])

        # fit
        guess = [I1_std, I1_amp, I1_mean-I1_std, I1_std, I1_amp]
        params, cov = sc.curve_fit(self._bimodal_fix_mean1, x, y, p0=guess)
        self.rI2mean = params[2]
        return params, cov


    '''
    Plotting functions
    '''
    def plot_calib_blobs(self):
        plt.figure()
        plt.plot(self.I1, self.Q1, '.', markersize='2', label=self.states[0])
        plt.plot(self.I2, self.Q2, '.', markersize='2', label=self.states[1])
        plt.axis('equal')
        plt.title('Calibration IQ blobs')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.tight_layout()

    def plot_rot_calib_blobs(self, plot_threshold=True):
        plt.figure()
        plt.plot(self.rI1, self.rQ1, '.', markersize='2', label=self.states[0])
        plt.plot(self.rI2, self.rQ2, '.', markersize='2', label=self.states[1])
        if plot_threshold:
            plt.axvline(self.threshold, color='k', linestyle='--')
        plt.axis('equal')
        plt.title('Rotated Calibration IQ blobs')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.tight_layout()

    def plot_calib_hist(self, bimodal1=True,bimodal2=True, plot_fit=True):
        plt.figure()
        plt.hist(self.rI1, bins=200, alpha=0.5, density=False, label=self.states[0])
        plt.hist(self.rI2, bins=200, alpha=0.5, density=False, label=self.states[1])

        if plot_fit:
            x1 = np.linspace(self.rI1.min(), self.rI1.max(), 200)
            x2 = np.linspace(self.rI2.min(), self.rI2.max(), 200)
            if bimodal1:
                P1, P2, params, cov = self.get_pop_from_fit(self.I1, self.Q1) # Do two-component fit but only use one component
                plt.plot(x1, params[2]*self._gauss(x1, params[0], params[1]), color='C0')
                self.rI1mean = params[0]  # set the mean based on only the correct component
            else:
                params, cov = self.get_gauss_fit(self.I1, self.Q1)
                self.rI1mean = params[0]
                plt.plot(x1,self._gauss(x1, params[0],params[1]), color='C0')#This will be incorrectly normalized but that doesn't matter!
            if bimodal2:
                #P1, P2, params, cov = self.get_pop_from_fit(self.I2, self.Q2) # Do two-component fit but only use one component
                params, cov = self.get_bimodal_fit_fix_mean1()
                plt.plot(x2, params[4]*self._gauss(x2, params[2], params[3]), color='C1')
                self.rI2mean = params[2] # set the mean based on only the correct component
            else:
                params, cov = self.get_gauss_fit(self.I2, self.Q2)
                self.rI2mean = params[0]
                plt.plot(x2,self._gauss(x2,params[0],params[1]), color='C1')#This will be incorrectly normalized but that doesn't matter!


        plt.xlabel('Rotated I')
        plt.legend()
        plt.tight_layout()

    def plot_blob(self, I, Q, label=None):
        plt.figure()
        plt.plot(I, Q, '.', markersize='2', label=label)
        plt.axis('equal')
        plt.title('IQ blobs')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.tight_layout()

    def plot_rot_blob(self, I, Q, label=None):
        rI, rQ = self._rot_blob(I, Q)
        plt.figure()
        plt.plot(rI, rQ, '.', markersize='2', label=label)
        plt.axis('equal')
        plt.title('Rotated IQ blobs')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.tight_layout()

    def plot_hist(self, I, Q, params, plot_fit=True,title=None):
        rI, rQ = self._rot_blob(I, Q)
        plt.figure()
        plt.hist(rI, bins=200, alpha=0.5, density=False)

        if plot_fit:
            x = np.linspace(rI.min(), rI.max(), 200)
            if len(params) == 6:
                p1_mean, p1_std, p1_amp, p2_mean, p2_std, p2_amp = params[0], params[1], params[2], params[3], params[4], params[5]
                P1 = (params[1] * params[2]) / (params[1] * params[2] + params[4] * params[5])
                P2 = (params[4] * params[5]) / (params[1] * params[2] + params[4] * params[5])
            else:
                p1_std, p1_amp, p2_std, p2_amp = params[0], params[1], params[2], params[3]
                p1_mean = self.rI1mean
                p2_mean = self.rI2mean
                P1 = (params[0] * params[1]) / (params[0] * params[1] + params[2] * params[3])
                P2 = (params[2] * params[3]) / (params[0] * params[1] + params[2] * params[3])

            if not math.isnan(P1):
                plt.plot(x, p1_amp * self._gauss(x, p1_mean, p1_std), linestyle='dotted', color='r', label=self.states[0])
                plt.plot(x, p2_amp * self._gauss(x, p2_mean, p2_std), linestyle='dashed', color='r', label=self.states[1])
                plt.plot(x, self._bimodal(x, p1_mean, p1_std, p1_amp, p2_mean, p2_std, p2_amp), linestyle='solid', color='r')

        if title is not None:
            plt.title(title)
        plt.xlabel('Rotated I')
        plt.legend()
        plt.tight_layout()


    def plot_hist_fix_means(self, I, Q, params, plot_fit=True):
        rI, rQ = self._rot_blob(I, Q)
        plt.figure()
        plt.hist(rI, bins=200, alpha=0.5, density=False)

        if plot_fit:
            x = np.linspace(rI.min(), rI.max(), 100)
            p1_std, p1_amp, p2_std, p2_amp = params[0], params[1], params[2], params[3]
            P1 = (params[0] * params[1]) / (params[0] * params[1] + params[2] * params[3])
            P2 = (params[2] * params[3]) / (params[0] * params[1] + params[2] * params[3])

            if not math.isnan(P1):
                plt.plot(x, 2 * p1_amp * self._gauss(x, np.mean(self.rI1), p1_std), linestyle='dotted', color='r', label='PL = '+str(round(P1, 3)))
                plt.plot(x, 2 * p2_amp * self._gauss(x, np.mean(self.rI2), p2_std), linestyle='dashed', color='r', label='PR = '+str(round(P2, 3)))
                plt.plot(x, 2 * self._bimodal_fix_means(x, p1_std, p1_amp, p2_std, p2_amp), linestyle='solid', color='r')

        plt.xlabel('Rotated I')
        plt.title('Fix Means')
        plt.legend()
        plt.tight_layout()