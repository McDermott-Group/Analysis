# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 13:44:03 2014

@author: Felix Jaeckel <felix.jaeckel@wisc.edu>
"""

import numpy as np
import matplotlib.pyplot as mpl

degree = np.pi/180.

class FourPointProbe():
    Recipe_4in_49pt_r = np.array([0, 13, 13, 13, 13, 13, 13, 13, 13, 26, 26., 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])
    Recipe_4in_49pt_theta = np.array([ 180, 180, 135,  90,  45,   0, 315, 270, 225, 180, 156, 135, 111,  90,  66,  45,  21,   0, 336, 315, 291, 270, 246, 225, 201, 180, 165, 149, 135, 120, 104,  90,  75,  59,  45,  30,  14,   0, 345, 329, 315, 300, 284, 270, 255, 239, 225, 210, 194])
    def __init__(self, recipe, fileName):
        if recipe == 'Felix 4in-49pt':
            d = np.loadtxt(fileName, delimiter=',')
            self.index = d[:,0] # "Point Num"
            self.r = d[:,1] # "mm fm Center"
            self.theta = d[:,2] # "Angle"
            self.rho = d[:,3] # "Ohms/Square"
            self.sym = d[:,4] # "ratio"
            self.currentRange = d[:,5] # "Current Range"
        elif recipe == 'Felix 4in-81pt':
            d = np.loadtxt(fileName, delimiter=',')
            self.index = d[:,0] # "Point Num"
            self.r = d[:,1] # "mm fm Center"
            if len(fileName) >= 14:
                ts = fileName[-14:]
                if ts < '201503130000000':
                    self.theta = d[:,2] # "Angle"
                    self.rho = d[:,3] # "Ohms/Square"
                    self.sym = np.zeros_like(self.rho) # Not in file
                    self.currentRange = d[:,4] # "Current Range"
                else:
                    self.theta = d[:,3] # "Angle"
                    self.rho = d[:,6] # "Ohms/Square"
                    self.sym = np.zeros_like(self.rho) # Not in file
                    self.currentRange = np.zeros_like(self.rho) # "Current Range"
        elif recipe == 'Felix Grid':
            d = np.loadtxt(fileName, delimiter=',')
            self.index = d[:,0] # "Point Num"
            self.r = d[:,1] # "mm fm Center"
            self.theta = d[:,2] # "Angle"
            self.rho = d[:,3] # "Ohms/Square"
            self.currentRange = d[:,4] # "Current Range"
            self.sym = d[:,5] # "ratio"
        elif recipe == '4in-49pt':
            d = np.loadtxt(fileName, delimiter=',')
            n = len(d)
            self.index = np.arange(1,n+1)
            self.r = self.Recipe_4in_49pt_r[:n] # Not available in the file, infer from similar program
            self.theta = self.Recipe_4in_49pt_theta[:n] # Not available in the file, infer from similar program
            self.rho = d
            self.sym = np.zeros_like(self.rho) # Not available in the file
            self.currentRange = -1 * np.ones_like(self.rho) # Not available in the file
        else:
            raise Exception("Unsupported recipe: %s" % recipe)

        if len(fileName) >= 14:
            ts = fileName[-14:]
        if ts < '20140715135335':
            # Before that time, the ROTATIONOFFSET=22.5 entry was missing in the AUTOEXEC.BAT file
            offset = 22.5 # According to AUTOEXEC.BAT in AFPP manual
        else: # Later files should be OK!
            offset = 0

        self.theta = 270.-offset-self.theta # Fix the default coordinate system

        self.theta *= degree

        self.x = np.cos(self.theta) * self.r
        self.y = np.sin(self.theta) * self.r

if __name__ == '__main__':

    Angstrom = 1E-10;
    cm = 1E-2;
    uOhm = 1E-6;

    from scipy.interpolate import griddata

    t = None
    rotate90deg = False

    wafer = '2696'; recipe = 'Felix 4in-49pt'; timestamp = '20170307115909'

    fpp = FourPointProbe(recipe, 'afppdata/%s/%s' % (recipe, timestamp))

    x = np.linspace(-50, +50, 501)
    y = np.linspace(-50, +50, 501)
    X,Y = np.meshgrid(x,y)


    rho = griddata((fpp.x, fpp.y), fpp.rho, (X, Y), method='cubic')

    extent=(-50,50, -50,50)

    # Construct wafer geometry
    onWafer = (X**2+Y**2) < 50.0**2
    onWafer &= Y > -45

    if rotate90deg:
        rho = np.fliplr(np.transpose(rho))
        onWafer = np.fliplr(np.transpose(onWafer))

    t=None
    if t is None:
        rhoUnits = u'$\Omega$/▫'
        rhoScale = 1
        thickness = ''
    else:
        rhoUnits = u'$\mu\Omega\cdot$cm'
        rhoScale = t/(uOhm*cm)
        thickness = '(t=%.1f$\AA$)' % (t/Angstrom)


    rho *= rhoScale
    fpp.rho *= rhoScale


    nPoints = len(fpp.rho)
    rhoMin = np.nanmin(fpp.rho)
    rhoMax = np.nanmax(fpp.rho)
    rhoMean = np.nanmean(fpp.rho)

    # At center
    r = np.sqrt(fpp.x**2+fpp.y**2)
    i = np.amin(r)
    rhoCenter = fpp.rho[i]
    print "At center:", rhoCenter, "Ohm/sq"
    print "Min:", rhoMin
    print "Mean:", rhoMean
    print "Max:", rhoMax
    #exit()

    #title = '%s %s- recipe: %s (%d points)\n' % (wafer, thickness, recipe, nPoints)
    #title += u'mean: %.3f %s (%.3f - %.3f %s)' % (rhoMean, rhoUnits, rhoMin, rhoMax, rhoUnits)
    title = u'Wafer %s\nmean %.3f %s (%.3f - %.3f %s)' % (wafer, rhoMean, rhoUnits, rhoMin, rhoMax, rhoUnits)
    mpl.figure()
    mpl.suptitle(title)
    mpl.imshow(onWafer, extent=extent, origin='lower')
    im = mpl.imshow(rho, extent=extent, origin='lower') # , vmin=0.131, vmax=0.136)
    CS = mpl.contour(X, Y, rho, colors='w', hold='on', origin='lower')
    mpl.scatter(fpp.x, fpp.y, s=2)
    mpl.clabel(CS, inline=1, fontsize=10, fmt = '%.2f')
    mpl.xlabel('x [mm]')
    mpl.ylabel('y [mm]')
    cb = mpl.colorbar(im)
    cb.set_label(rhoUnits)
    mpl.savefig('%s_Abs.pdf' % wafer)
    mpl.savefig('%s_Abs.png' % wafer)

    normalized = 100. * rho / np.nanmin(rho)

    mpl.figure()
    mpl.suptitle(title)
    mpl.imshow(onWafer, extent=extent, origin='lower')
    im = mpl.imshow(normalized, extent=extent, origin='lower')
    CS = mpl.contour(X, Y, normalized, colors='w', hold='on', origin='lower')
    mpl.scatter(fpp.x, fpp.y, s=2)
    mpl.clabel(CS, inline=1, fontsize=10, fmt = '%.0f')
    mpl.xlabel('x [mm]')
    mpl.ylabel('y [mm]')
    cb = mpl.colorbar(im)
    cb.set_label(u'%')
    mpl.savefig('%s_Percent.pdf' % wafer)
    mpl.savefig('%s_Percent.png' % wafer)
    mpl.show()
