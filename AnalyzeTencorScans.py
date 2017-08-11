# -*- coding: utf-8 -*-
"""
Plot Tencor scan data and extract curvature, stress, etc.
Created on Fri Mar 14 13:22:33 2014
@author: Felix Jaeckel <jaeckel2@wisc.edu>
"""

import matplotlib.pyplot as mpl
import numpy as np
from Tencor import ScanImporter


def labelForScan(scan, i=None):
    if i is not None:
        label = 'scan %i %s' % (i+1, scan.comment) # Have to add an extra one here because humans count from 1
    else: # Not sure why one would want sample ID here, maybe if that was used to distinguish between before and after scans?
        label='%s %s' % (scan.sampleId, scan.comment)
    return label


def remap(x1, y1, x2):
    if np.allclose(x1,x2):
        return y1
    else:
        return np.interp()

def commonSpan():
    pass

def stoneyStress(Ebiaxial, tWafer, tFilm, DeltaR):
    stress = (tWafer**2 / tFilm) * (1./(6*DeltaR)) * Ebiaxial
    return stress

def analyzeScans(scans, minX=-1E10, maxX=+1E10, minIntensity=0, wafer=''):
    nScans = len(scans)
    sampleId = scans[0].sampleId
    print "Sample ID:", sampleId
    assert(nScans > 0)
    waferThickness = scans[0].waferThickness
    #waferThickness = 3E-4  # Hack!!!
    print "Wafer thickness:", waferThickness / 1E-6, 'um'
    modulus = scans[0].modulus
    print "Modulus:", modulus

    for scan in scans:
        if scan.waferThickness != waferThickness:
            print 'Wafer thickness does not match!'
        if scan.modulus != modulus:
            print 'Wafer modulus does not match!'

    colors = ['r', 'g', 'b', 'y', 'm']

    nPlots = 2
    if nScans > 1:
        nPlots += nScans-1
    nPlot = 1
    mpl.figure()
    mpl.subplot(nPlots, 1, nPlot); nPlot += 1
    mpl.title('Wafer: %s/%s' % (wafer, sampleId))
    mpl.ylabel('Intensity [V]')
    #mpl.xlabel('x [mm]')

    for i,scan in enumerate(scans):
        mpl.plot(1E3*scan.x, scan.intensity, '.-', color=colors[i], label=labelForScan(scan,i))

    legend = mpl.legend(loc='lower right')
    legend.get_frame().set_alpha(0.5)

    mpl.subplot(nPlots, 1, nPlot); nPlot += 1
    mpl.ylabel('Deflection [rad]')
   # mpl.xlabel('x [mm]')

    for i, scan in enumerate(scans):
        print 'Film', i, 'Bow:', scan.maxBow / 1E-6, 'um'
        iFit = (scan.x >= minX) & (scan.x <= maxX) & (scan.intensity > minIntensity)
        fit,pcov = np.polyfit(scan.x[iFit], scan.deflection[iFit], 1, cov=True)
        slope = fit[0]
        slopeStd = np.sqrt(np.diag(pcov))[0]
        R = 1./slope
        RStd = R*(slopeStd/slope)
        mpl.plot(1E3*scan.x, scan.deflection, '.', color=colors[i])
        mpl.plot(1E3*scan.x[iFit], np.polyval(fit, scan.x[iFit]), '-', color = colors[i], label='$R=%.2f\pm%.2f m (%.2f m)$' % (R, RStd, scan.radius))
    legend = mpl.legend(loc='lower right')
    legend.get_frame().set_alpha(0.5)

    # Now analyze in pairs
    stresses = []
    for i in range(nScans-1):
        scanA = scans[i]
        scanB = scans[i+1]
        mpl.subplot(nPlots, 1, nPlot); nPlot += 1
        mpl.ylabel('Bow [um]')
        mpl.plot(1E3*scanA.x, 1E6*scanA.bow, '.', color = colors[i], label=labelForScan(scan,i))
        mpl.plot(1E3*scanB.x, 1E6*scanB.bow, '.', color = colors[i+1], label=labelForScan(scan,i+1))

        assert(np.allclose(scanA.x, scanB.x))
        x = scanA.x
        y = scanB.bow - scanA.bow
        intensity = np.min(np.vstack([scanA.intensity, scanB.intensity]), axis=0)
        iFit = (x >= minX) & (x <= maxX) & (intensity > minIntensity)
        fit,pcov = np.polyfit(x[iFit], y[iFit], 2, cov=True)
        a = fit[0]
        aStd = np.sqrt(np.diag(pcov))[0]
        DeltaR=1./(2*a)
        DeltaRStd = R*(aStd/a)
        #scanB.filmThickness=400E-10
        print "Substrate modulus:", modulus
        if scanB.filmThickness == 0:
            print "### Warning: film thickness was 0!"
            scanB.filmThickness = 7654.0

        print "Film %i thickness: %.1f A" % (i, scanB.filmThickness / 1E-10)
        stress = stoneyStress(modulus, waferThickness, scanB.filmThickness, DeltaR)
        stressStd = stress*(aStd/a)
        print 'Film %i stress: %.2f MPa' % (i, 1E-6*stress)
        mpl.plot(1E3*x, 1E6*y, '.', color='k', label = 'scan %i - %i' % (i+2, i+1)) # Have to add an extra one here because humans count from 1
        mpl.plot(1E3*x[iFit], 1E6*np.polyval(fit,x[iFit]), '-', color='k', label='$\Delta R=%.2f\pm%.2f m$' % (DeltaR, DeltaRStd) )
        legend = mpl.legend(loc='lower right')
        legend.get_frame().set_alpha(0.5)
        mpl.title('Film #%d (%g A): Stress $%g\pm%g$ MPa' % (i+1, 1E10*scanB.filmThickness, 1E-6*stress, 1E-6*stressStd))
        stresses.append(stress)

    mpl.xlabel('x [mm]')
    mpl.show()
    return stresses


if __name__ == '__main__':
#    import glob
    mpl.close('all')

    case = 4
    if case == 1:
        fileName = 'KLK/2501.SCN'
        scans = ScanImporter(fileName)
        print "Number of scans in file:", scans.numberOfScans()
        minX = 30E-3
        maxX = 70E-3
        #analyzeScans(scans, minX, maxX)
        analyzeScans(scans)
        #scans = [scans[3], scans[4]]  # For 2402_30, 60
    elif case == 2:
        scans = ScanImporter('FJ2/FJ2.SCN')
        scans = scans['1013']
    elif case == 3:
        scans = ScanImporter('KLK2/Stress20140710/2422_180.SCN')
        print "Number of scans in file:", scans.numberOfScans()
    elif case == 4: # When all the scans are together in one file
        wafer = '2699'
        fileName = '/Volumes/smb/mcdermott-group/data/filmStress/KLK_testData/%s.SCN' % wafer
        scans = ScanImporter(fileName)
        print "Number of scans in file:", scans.numberOfScans()
        sampleIds = scans.sampleIds
        uniqueIds = set(sampleIds)
        angles = np.asarray([float(uid) for uid in uniqueIds])
        filmStress = []
        for uid in uniqueIds:
            iForAngle = np.where(np.array(sampleIds) == uid)[0]
            scansForAngle = []
            for i in iForAngle:
                scan = scans[int(i)]
                print "Scan angle:", scan.angle
                scansForAngle.append(scan)
            if True: #'27' in uid:
                #scansForAngle = [scansForAngle[0], scansForAngle[1], scansForAngle[3]]
                stresses = analyzeScans(scansForAngle, wafer=wafer) #[1:3])
                filmStress.append(stresses[0]) # 0 or 1 depending which layer
                #print filmStress
            #d = raw_input('Enter')
            #mpl.close('all')

        filmStress = np.asarray(filmStress)
        iSort = np.argsort(angles)
        print "Angles:", angles[iSort]
        print "Film stress (MPa):", '\t'.join(['%.3f' % x for x in filmStress[iSort]/1E6])

