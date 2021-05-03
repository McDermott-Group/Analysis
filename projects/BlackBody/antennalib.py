"""
A library for the antenna model
"""

import numpy as np
from scipy.constants import *

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