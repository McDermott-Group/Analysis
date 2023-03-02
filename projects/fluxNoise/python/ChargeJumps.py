import numpy as np
import noiselib
import importlib
importlib.reload(noiselib)
    
class ve(object):

    def __init__(self, v, e):
        self.v = v
        self.e = e
    
    def __str__(self):
        sig_dig = int(np.floor(np.log10(np.abs(self.e))))
        return '{:.{sd}f}({})'.format( self.v, 
                                       int(np.round(10**(-sig_dig) * self.e)), 
                                       sd=-sig_dig )
    
    def __repr__(self):
        return str(self.to_tuple(self))
    
    def __float__(self):
        return float(self.v)
    
    def __int__(self):
        return int(self.v)
    
    def __abs__(self):
        return ve(np.abs(self.v),self.e)
    
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
        q = noiselib.alias(q, 0.5)
        return 1. * pve(np.sum(np.abs(q)>thresh)) / pve(q.size)
    
    def g_thresh(self, q, thresh, dt):
        q, dt = noiselib.overlap(q, dt)
        q = noiselib.alias(q, 0.5)
        return 1. * pve(np.sum(np.abs(q)>thresh)) / np.sum(dt)
    
    def thresh_fraction(self, q, thresh):
        q = q[np.isfinite(q)]
        q = noiselib.alias(q, 0.5)
        return 1. * pve(np.sum(np.abs(q)>thresh)) / pve(q.size)
    
    def assym(self, q, thresh):
        """ """
        q = q[np.isfinite(q)]
        q = noiselib.alias(q, 0.5)
        p = pve( np.sum( q > thresh ) )
        a = pve( np.sum( np.abs(q) > thresh ) )
        return 1. * p / a
    
    def assym1324(self, q1, q2, thresh1, thresh2):
        """q1 and q2 should be arrays of charge measurements, with equal lengths.
        nans can be included and those rows will be removed from both arrays."""
        q1, q2 = noiselib.overlap(q1, q2)
        q1 = noiselib.alias(q1, 0.5)
        q2 = noiselib.alias(q2, 0.5)
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
        q1 = noiselib.alias(q1, 0.5)
        q2 = noiselib.alias(q2, 0.5)
        n_corr = pve( np.sum( ( np.abs(q1) > thresh1 ) & ( np.abs(q2) > thresh2 ) ) )
        print(1.*n_corr/pve( np.sum( ( np.abs(q1) > thresh1 ) ) ))
        print(1.*n_corr/pve( np.sum( ( np.abs(q2) > thresh2 ) ) ))
        pABo = 1.*n_corr / pve(q1.size)
        return pABo
    
    def raw_correlation(self, q1,q2, thresh1, thresh2):
        """The number of corrrelated events divided by the number of events in 
        which either has a jump above the threshhold."""
        corr_counts = pve(np.sum( (np.abs(q1)>thresh1) & (np.abs(q2)>thresh2) ))
        tot_counts = pve(np.sum( (np.abs(q1)>thresh1) | (np.abs(q2)>thresh2) ))
        return 1. * corr_counts / tot_counts
        
    def calc_probabilities(self, qA, qB, threshA, threshB):
        
        pAo = self.p_thresh(qA, threshA)
        pBo = self.p_thresh(qB, threshB)
        pABo = self.p_correlation(qA, qB, threshA, threshB)
        
        pAB = 1.*(pABo-pAo*pBo) / (1.+pABo-pAo-pBo)
        pA = (pAo-pAB) / (1.-pAB)
        pB = (pBo-pAB) / (1.-pAB)
        pAB_norm = 1.*pAB / ((pAo+pBo)/2.)
        
        return pAo, pBo, pA, pB, pABo, pAB, pAB_norm
        
    def calc_rates(self, qA, qB, threshA, threshB, dt):
        
        _,_,_,_,_,_,pAB_norm = self.calc_probabilities(qA, qB, threshA, threshB)
        gAo = self.g_thresh(qA, threshA, dt)
        gBo = self.g_thresh(qB, threshB, dt)
        gABo = (gAo+gBo)/2.
        gAB = pAB_norm * gABo
        
        return gAo, gBo, gAB
        
    def calc_params(self, qA, qB, threshA, threshB, dt):
    
        a1324 = self.assym1324(qA, qB, threshA, threshB)
        pAo, pBo, pA, pB, pABo, pAB, pAB_norm = \
            self.calc_probabilities(qA, qB, threshA, threshB)
        gAo, gBo, gAB = self.calc_rates(qA, qB, threshA, threshB, dt)
        # return pAB_norm, a1324, gAo, gBo
        return pAo, pBo, pA, pB, gAo, gBo, pABo, pAB, gAB, pAB_norm, a1324
        
