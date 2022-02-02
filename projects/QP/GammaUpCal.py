from antennalib import getGammaUpDownFromQPSpectrum, XqpUpDown
import matplotlib.pyplot as plt
import numpy as np

if 1:
    QPSpectrumFile = "QPEnergySpectrum_2msec.txt"
    # QPSpectrumFile = "QPEnergySpectrum_20usec.txt"
    QPSpectrum = np.loadtxt(QPSpectrumFile)
    # print('QPSpectrum=', QPSpectrum[:,0])
    epsilon = QPSpectrum[:, 0]
    n = QPSpectrum[:, 1]

    gammaup = getGammaUpDownFromQPSpectrum(epsilon, n)[0]
    gammadown_ = getGammaUpDownFromQPSpectrum(epsilon, n)[1]

    gammadown = (1.0/(1.414*1e4*5e-15)) * 10 ** 1.5 * 2 * sum(n)

    print ('gammaup=', gammaup/1e6, 'MHz')
    print ('gammadown_=', gammadown_/1e6, 'MHz')
    print ('gammadown=', gammadown/1e6, 'MHz')
    print('x_qp=', 2*sum(n))

if 0:
    xQP = XqpUpDown()
    T1 = [50.0e-6, 5.0e-6, 10.0e-6]
    xQP.import_data(T1)