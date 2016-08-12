import scipy as sp
from scipy import *
from scipy.optimize import curve_fit
import numpy as np
from numpy import *
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

########################################################################################################################
########################################################################################################################
###     The following code reads in .S2P files from "saveFitResultsPath". The fitting expects the data to be
###     freq, real, and imaginary. The columns to be read in are specified in "dataColumns." The function S21Fit
###     is taken from Ben Mazin's thesis.
###
###
###
###


# Path to save the fit results to
saveFitResultsPath = '/Users/matthewbeck/Box Sync/Data/multiplexed resonators/Dip Tests/10-15-2015/30 um Trace/19-24 ums/'

# fit results file name
fitResultsFileName = 'fits.dat'

# text header for fit results file
textHeader = "f0 \t  df0\t\tQi \t      dQi \t    Qc \t      dQc"

# path to save fit figures
saveFigurePath = '/Users/matthewbeck/Box Sync/Data/multiplexed resonators/Dip Tests/10-15-2015/30 um Trace/19-24 ums/'

# fit figures file name
outputFigureName = "Fits.pdf"

# initiliaze output figure to write to
outputFigure = PdfPages(saveFigurePath + outputFigureName)

# Data path
dataPath = '/Users/matthewbeck/Box Sync/Data/multiplexed resonators/Dip Tests/10-15-2015/30 um Trace/19-24 ums/'

# data extension
dataExt = '*.s2p'

def S21Fit(f, f0, Qi, Qc ,L, through, alpha):
    # minimum value at the dip
    sMin = Qc/(Qi + Qc)

    # Total quality factor
    qTotal = Qi*Qc/(Qi + Qc)

    # Normalized frequency difference
    dx = (f - f0)/f0

    # First data point in frequency array
    ff = f[0]

    # Mazin S21 function, the added 1i*L term is for dip assymetry about f0
    s21 = (sMin + 2*1j*qTotal*dx)/(1 + 2*1j*qTotal*dx + 1j*L)

    # logMag of S21 function that also accounts for through loss and angle (alpha) of the curve
    s21LogMag = (alpha*(ff-f))+10 * sp.log10(s21*np.conj(s21)) + through
    return s21LogMag.real

# Data Path
path = dataPath + dataExt

# number of header rows in S2P file
headerRows = 0

# Columns from which to read in data
dataColumns = (0,1,2)

# list of file names in "path"
files = glob.glob(path)

# initialize i -> file counter
i = 0

# initialize data outputArray
outputArray = []

# lower Bounds on fitting params
QiMin = 1e3
QcMin = 1e3
lMin = -1
alphaMin = -.5

# upper Bounds on Fitting params
QiMax = 1e5
QcMax = 1e5
lMax = 1
alphaMax = .5

# While loop to go through files and fit them
while i < len(files):

    # write in data from file
    f,realData,imagData = sp.loadtxt(str(files[i]),skiprows=headerRows, usecols = dataColumns, unpack=True)

    # combine real and complex data
    complexData = realData + 1j*imagData

    # logMag of complexData
    logMagData = (10*sp.log10(complexData*sp.conj(complexData))).real

    # Find minimum value of logMagData and its array index
    val, idx = min((val, idx) for (idx, val) in enumerate(logMagData))

    # logMagData thru value estimation
    thru = sp.average(logMagData[0:10])

    # bounds for fitting algorithm
    #            [fmin, QiMin, QcMin, Lmin, thruMin, alphaMin]
    lowerBound = [f[idx]-5e6, QiMin, QcMin, lMin, thru-5, alphaMin]

    #            [fMax, QiMax, QcMax, lMax, thruMax, alphaMax]
    upperBound = [f[idx]+5e6, QiMax, QcMax, lMax, thru+5, alphaMax]

    # Invoke fitting routine
    # popt are the optimized fitting parameters.
    # pcov is the convariance matrix
    popt, pcov = curve_fit(S21Fit,f,logMagData,bounds=(lowerBound, upperBound))

    # 1 sigma error is squareroot of diagonal of covariance matrix
    perr = np.sqrt(np.diag(pcov))


    #fits = sp.asarray([popt, perr]).T

    # Write fits and error to outputArray
    outputArray.append(
        [
            popt[0], # fFit
            perr[0], # dF
            popt[1], # QiFit
            perr[1], # dQi
            popt[2], # QcFit
            perr[2]  # dQc
        ]
    )

    # Calculate S21 profile with fitted parameters
    fitOutPut=S21Fit(f,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

    # Plotting
    fig, ax = plt.subplots()

    # plot Data, x-axis is in units of GHz
    plt.plot(f/1e9, logMagData, 'bo', label = "Data")

    # plot fit + Qi, Qc #'s & error
    plt.plot(f/1e9, fitOutPut, 'r--',linewidth = 2, label = "Fit \n" +
            "Qi = " + str(round(popt[1],0)) + " +/- " + str(round(perr[1],0)) + "\n"+
            "Qc = " + str(round(popt[2],0)) + " +/- " + str(round(perr[2],0))
             )

    # put legend on graph
    plt.legend()

    # set y axis label
    plt.ylabel('S21')

    # set x axis label
    plt.xlabel('Frequency (GHz)')

    # turn grid lines on
    plt.grid(b=True, which='major', color='gray', linestyle='--')

    # put legend in bottom left hand corner
    ax.legend(loc='lower left')

    # strech axis back a little
    xAxisMin = f[0]/1e9 - .005
    xAxisMax = f[-1]/1e9

    # Change limits of x axis so legen properly fits
    plt.xlim(xAxisMin,xAxisMax)

    # save fit figure to outputFigure
    outputFigure.savefig()

    # update counter
    i = i + 1

# Close the output figure
outputFigure.close()

# write out the fit data
np.savetxt(saveFitResultsPath+fitResultsFileName, outputArray,fmt='%.3e', header = textHeader)