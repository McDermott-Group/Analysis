import numpy as np
import matplotlib.pyplot as plt

class CST_Impedance:
    f = 0
    z_real = 0
    z_img = 0

    def __init__(self, filename, frequency_list=[]):
        tmp_f = np.loadtxt(filename, usecols=[0], skiprows=3)
        tmp_real = np.loadtxt(filename, usecols=[1], skiprows=3)
        tmp_img = np.loadtxt(filename, usecols=[2], skiprows=3)

        if len(frequency_list) > 0:
            self.z_real = np.interp(frequency_list,tmp_f,tmp_real)
            self.z_img = np.interp(frequency_list,tmp_f,tmp_img)
            self.f = 1e9 * tmp_f
        else:
            self.z_real = tmp_real
            self.z_img = tmp_img
            self.f = 1e9 * tmp_f

class Antenna:

    r = 0
    c = 0
    z_real = 0
    z_img = 0

    def __init__(self, R, C, z_real, z_img,):
        pass

    def e_c(self):
        pass

class Radiator(Antenna):

    pass