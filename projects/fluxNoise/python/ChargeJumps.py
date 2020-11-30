import os
import numpy as np
import noiselib
reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from scipy.stats import norm
from dateStamp import dateStamp
import glob

# def add(*args):
    # res = np.sum([v for v,e in args])
    # err = np.sqrt(np.sum([e**2 for v,e in args]))
    # return res, err

# def mult(*args):
    # res = np.product([v for v,e in args])
    # err = np.abs(res) * np.sqrt(np.sum([(e/v)**2 for v,e in args]))
    # return res, err

# def div((A,dA),(B,dB)):
    # return mult( (A,dA), (1./B,dB/B**2) )
    
class ve(object):

    def __init__(self, v, e):
        self.v = v
        self.e = e
    
    def __str__(self):
        print 'this was called'
        return str(self.__repr__())
    
    def __repr__(self):
        return str(self.to_tuple(self))
    
    def to_tuple(self, o):
        if isinstance(o, ve):
            return (o.v, o.e)
        elif isinstance(o, (tuple,list)):
            return tuple(o)
        else:
            return (o,0.0)

    def __add__(self, o):
        v,e = self.to_tuple(o)
        res = self.v + v
        err = np.sqrt(self.e**2 + e**2)
        return ve(res, err)
    __radd__ = __add__

    def __sub__(self, o):
        v,e = self.to_tuple(o)
        res = self.v - v
        err = np.sqrt(self.e**2 + e**2)
        return ve(res, err)
    def __rsub__(self, o):
        v,e = self.to_tuple(o)
        res = v - self.v
        err = np.sqrt(self.e**2 + e**2)
        return ve(res, err)

    def __mul__(self, o):
        v,e = self.to_tuple(o)
        res = self.v * v
        err = np.abs(res) * np.sqrt((self.e/self.v)**2 + (e/v)**2)
        return ve(res, err)
    __rmul__ = __mul__

    def __truediv__(self, o):
        try:
            v,e = self.to_tuple(o)
            res = 1. * self.v / v
            err = np.abs(res) * np.sqrt((self.e/self.v)**2 + (e/v)**2)
            return ve(res, err)
        except ZeroDivisionError:
            return ve(np.nan, np.nan)
    __div__ = __truediv__
    def __rtruediv__(self, o):
        try:
            v,e = self.to_tuple(o)
            res = 1. * v / self.v
            err = np.abs(res) * np.sqrt((self.e/self.v)**2 + (e/v)**2)
            return ve(res, err)
        except ZeroDivisionError:
            return ve(np.nan, np.nan)
    __rdiv__ = __rtruediv__

    def __pow__(self, o):
        v,e = self.to_tuple(o)
        res = self.v ** v
        err = np.abs( res*v*self.e/self.v )
        return ve(res, err)
    def __rpow__(self, o):
        v,e = self.to_tuple(o)
        res = v ** self.v
        err = np.abs( res*self.v*e/v )
        return ve(res, err)
    

def pve(v):
    """For Poisson processes, the error in number of counts is sqrt(counts)."""
    return ve(v,np.sqrt(v))
    

class ChargeJumps(object):
    
    def set_charge(self, q_induced):
        
        self.q_induced = q_induced
    
    def p_thresh(self, q, thresh):
        q = q[np.isfinite(q)]
        return 1. * pve(np.sum(np.abs(q)>thresh)) / pve(q.size)
    
    def g_thresh(self, q, thresh, dt):
        q, dt = noiselib.overlap(q, dt)
        return 1. * pve(np.sum(np.abs(q)>thresh)) / np.sum(dt)
    
    def assym1324(self, q1, q2, thresh1, thresh2):
        """q1 and q2 should be arrays of charge measurements, with equal lengths.
        nans can be included and those rows will be removed from both arrays."""
        q1, q2 = noiselib.overlap(q1, q2)
        Q1 = pve( np.sum( (  q1 > thresh1 ) & (  q2 > thresh2 ) ) )
        Q2 = pve( np.sum( ( -q1 > thresh1 ) & (  q2 > thresh2 ) ) )
        Q3 = pve( np.sum( ( -q1 > thresh1 ) & ( -q2 > thresh2 ) ) )
        Q4 = pve( np.sum( (  q1 > thresh1 ) & ( -q2 > thresh2 ) ) )
        a1324 = 1.*(Q1+Q3)/(Q2+Q4)
        return a1324
    
    def p_correlation(self, q1, q2, thresh1, thresh2):
        """q1 and q2 should be arrays of charge measurements, with equal lengths.
        nans can be included and those rows will be removed from both arrays."""
        q1, q2 = noiselib.overlap(q1, q2)
        n_corr = pve( np.sum( ( np.abs(q1) > thresh1 ) & ( np.abs(q2) > thresh2 ) ) )
        p_obs = 1.*n_corr / pve(q1.size)
        return p_obs
        
        
    def plot_charge_correlation(self, CO, label1, label2, thresh=None, 
                                      datasets=None, ax=None, plot=True):
        
        jumps, sigma = CO.get_jump_sizes(datasets)
        dt = CO.get_time_steps(datasets)
        qA, qB = jumps[label1], jumps[label2]
        threshA, threshB = thresh
        
        pA = self.p_thresh(qA, threshA)
        pB = self.p_thresh(qB, threshB)
        gA = self.g_thresh(qA, threshA, dt)
        gB = self.g_thresh(qB, threshB, dt)
        gAB = (gA+gB)/2.
        
        a1324 = self.assym1324(qA, qB, threshA, threshB)
        p_obs = self.p_correlation(qA, qB, threshA, threshB)
        
        pC = 1.*(p_obs-pA*pB) / (1.+p_obs-pA-pB)
        pAp = (pA-pC)/(1.-pC)
        pBp = (pB-pC)/(1.-pC)
        pCp = 1.*pC/((pA+pB)/2.)
        
        gC = pCp * gAB
        
        return pCp, a1324, gA, gB
