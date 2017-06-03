# -*- coding: utf-8 -*-
"""
Read and digest Tencor scan (.SCN), data (.DAT), graph (.GRP), and deflection (.TXT) files.
Reverse engineered based on Tencor WinLFX 4.40 (registered trademark) software.
Although the binary file format is mostly decoded, the interpretation of some of the numbers still remains unclear.

@author: Felix Jaeckel <jaeckel2@wisc.edu>
Created on Thu Mar 13 10:21:19 2014

Tencor SCN (scan) file format (holds up to 1000 records):
Layout:
0x000 - 0x003 Magic header A5A5 C800
0x004 - 0x063 Mystery stuff, looks like uninitialized (sometimes see FAT16 records here)
0x064 - 0x256 First record, length 0x1F3
0x257 - 0x449 Second record, same length
0x44A -       Third record, etc...

Each record:
0x000 - 0x010  1*string Sample ID, 0x0 terminated (16 maximum + 0x0)
0x011 - 0x01D  1*string Comment, 0x0 terminated (12 maximum + 0x0
0x01E - 0x01F  1*short  Temperature [C]
0x020 - 0x023  1*float  Radius of curvature [m]
0x024 - 0x025  1*short  Film thickness [A]
0x026 - 0x027  2*bytes  MysteryS (00 00)
0x028 - 0x029  1*short  Wafer thickness [um]
0x02A - 0x02D  1*float  Elastic modulus [1E11 Pa]
0x02E - 0x02F  1*short  Diameter [mm]
0x030 - 0x031  1*short  Angle [degree]
0x032 - 0x035  1*float  Scan start [mm]
0x036 - 0x039  1*float  Scan end [mm]
0x03A - 0x03D  1*float  Scan step [cm] (Doesn't quite make sense!)
0x03E - 0x041  1*uint   Date/Time, seconds since 00:00:00 Jan 1st, 1970 (Unix epoch)
0x042 - 0x042  1*byte   First  (1=first, 0=single)
0x043 - 0x046  1*float  Working distance = 14.389 (Wafer/detector distance (14.41") - Wafer thickness (in inches))
0x047 - 0x04A  1*float  Bow [um]
0x04B - 0x04E  1*float  MysteryF2 = 20.367, also 22.918, 15.265, 15.897
0x04F - 0x052  1*float  Low intensity alarm [V]
0x053 - 0x11A 50*floats Detector signal A-B [V]
0x11B - 0x1E2 50*floats Detector signal A+B [V]
0x1E3 - 0x1E6  1*float  Mystery F3 E31B C8BA -0.00152671  Or together as...
0x1E7 - 0x1EA  1*float  Mystery F4 E72F 1B3F +0.60619968  ... double: 1.037E-4
0x1EB - 0x1EF all zeros

Not found:
- Laser (wavelength/mode)
- Stress
- Hole

DAT File:
No header, all records just one after another
Length of each block: 0x64 = 100 bytes
0x000 - 0x010  1*string Sample ID, 0x0 terminated (16 maximum + 0x0)
0x011 - 0x01D  1*string Comment, 0x0 terminated (12 maximum + 0x0
0x01E - 0x01F  1*short  Temperature [C]
0x020 - 0x023  1*float  Radius of curvature [m]
0x024 - 0x027  1*float  Stress [Pa]
0x028 - 0x029  1*short  Film thickness [A]
0x02A - 0x02B  1*short  00 00 ?
0x02C - 0x02D  1*short  Wafer thickness [um]
0x02E - 0x031  1*float  Elastic modulus [1E11 MPa]
0x032 - 0x033  1*short  Wafer diameter [mm]
0x034 - 0x037  1*float  Position (of what?!) [?]
0x038 - 0x03B  1*float  Intensity [V] (mean/max?)
0x03C - 0x03D  1*short  Angle [degree]
0x03E - 0x041  1*float  0
0x042 - 0x045  1*float  Scan start [mm]
0x046 - 0x049  1*float  Scan stop [mm]
0x04A - 0x04D  1*float  Scan step [cm] (Doesn't quite make sense!)
0x04E - 0x051  1*float  Bow [um]
0x052 - 0x055  1*uint   Date/Time, seconds since 00:00:00 Jan 1st, 1970 (Unix epoch)
0x056 - 0x057  1*short  00 02 flag? (sometimes 01 02)
0x058 - 0x05B  1*float  00 00 FF FF = nan always?
0x05C - 0x05F  1*float  Radius of curvature [m] (again!)
0x060 - 0x063  1*float  Bow [um] (again!)
"""


import numpy as np
import struct
import time
import scipy.integrate as integrate
import pandas as pd

class TencorTable(object):
    def __init__(self, fileName):
        df = pd.DataFrame(columns=['sampleId', 'comment', 'T','R','stress','film_t', 'wafer_t', 'wafer_E', 'wafer_d', 'mystery1', 'intensity', 'angle', 'mystery2', 'scanStart', 'scanStop', 'scanStep', 'bow', 'time', 'mystery3', 'mystery4'])
        with open(fileName, 'rb') as f:
            while True:
                record = f.read(100)
                #print len(record)
                if len(record) < 100:
                    break
                sampleId  = struct.unpack('17s', record[0x00:0x11])[0].split('\0', 1)[0]
                comment   = struct.unpack('13s', record[0x11:0x1E])[0].split('\0', 1)[0]
                T         = struct.unpack('<H', record[0x1E:0x20])[0]
                R         = struct.unpack('<f', record[0x20:0x24])[0]
                stress    = struct.unpack('<f', record[0x24:0x28])[0]
                film_t    = struct.unpack('<I', record[0x28:0x2C])[0] * 1E-1
#                mystery1  = struct.unpack('<H', record[0x2A:0x2C])[0]
                wafer_t   = struct.unpack('<H', record[0x2C:0x2E])[0] * 1E-6
                wafer_E   = struct.unpack('<f', record[0x2E:0x32])[0] * 1E11*1E6
                wafer_d   = struct.unpack('<H', record[0x32:0x34])[0] * 1E-3
                mystery1  = struct.unpack('<f', record[0x34:0x38])[0]
                intensity = struct.unpack('<f', record[0x38:0x3C])[0]
                angle     = struct.unpack('<H', record[0x3C:0x3E])[0]
                mystery2  = struct.unpack('<f', record[0x3E:0x42])[0]
                scanStart = struct.unpack('<f', record[0x42:0x46])[0] * 1E-3
                scanStop  = struct.unpack('<f', record[0x46:0x4A])[0] * 1E-3
                scanStep  = struct.unpack('<f', record[0x4A:0x4E])[0] * 1E-2
                bow       = struct.unpack('<f', record[0x4E:0x52])[0] * 1E-6
                time      = struct.unpack('<I', record[0x52:0x56])[0]
                mystery3  = struct.unpack('<H', record[0x56:0x58])[0]
                mystery4  = struct.unpack('<f', record[0x58:0x5C])[0]
                R2        = struct.unpack('<f', record[0x5C:0x60])[0]
                bow2      = struct.unpack('<f', record[0x60:0x64])[0] * 1E-6
                assert R == R2, 'The two radii of curvature do not match'
                assert bow == bow2, 'The two wafer bow records do not match'

                df.loc[len(df)] = {'sampleId': sampleId, 'comment':comment, 'T':T,'R':R,'stress':stress,'film_t':film_t,
                             'wafer_t':wafer_t, 'wafer_E':wafer_E, 'wafer_d':wafer_d,
                             'mystery1':mystery2, 'intensity':intensity, 'angle':angle,
                             'mystery2':mystery3, 'scanStart':scanStart, 'scanStop':scanStop, 'scanStep':scanStep,
                             'bow':bow, 'time':time, 'mystery4':mystery3, 'mystery5':mystery4}
        self.df = df

    def saveToExcel(self, fileName):
        t.df.to_excel(fileName)



def linearFit(x, y):
    return np.sum(x*y)/np.sum(x*x)

def loadGrpFile(fileName):
    '''Load data from a Tencor grp file.
    This is pretty generic in that it just ignores everything that does not conform to two numeric columns.'''
    with open(fileName, 'r') as f:
        lines = f.readlines()
        xs = []
        ys = []
        for l in lines:
            a = l.split()
            if len(a) != 2:
                continue
            try:
                x = float(a[0])
                y = float(a[1])
                xs.append(x)
                ys.append(y)
            except:
                pass
    return np.array(xs), np.array(ys)

class Struct:
    def __init__ (self, *argv, **argd):
        if len(argd):
            # Update by dictionary
            self.__dict__.update (argd)
        else:
            # Update by position
            attrs = filter (lambda x: x[0:2] != "__", dir(self))
            for n in range(len(argv)):
                setattr(self, attrs[n], argv[n])

class Wave(Struct):
    x = []
    y = []

def loadTripleGrpFile(fileName):
    x, y = loadGrpFile(fileName)
    assert(len(x)==len(y))
    assert(len(x) % 3 == 0)
    n = int(len(x) / 3)
    curve1 = Wave(x[0*n:1*n], y[0*n:1*n])
    curve2 = Wave(x[1*n:2*n], y[1*n:2*n])
    curve3 = Wave(x[2*n:3*n], y[2*n:3*n])
    return curve1, curve2, curve3

class Scan(object):
    '''Represents a Tencor stress measurement scan record.'''
    RECORD_LENGTH=0x1F3
    OFFSET_COMMENT=0x11
    OFFSET_TEMPERATURE=0x1E # short
    OFFSET_RADIUS = 0x20 # float?
    OFFSET_FILMTHICKNESS=0x24 # short
    OFFSET_MYSTERY_S = 0x26 # short?
    OFFSET_WAFERTHICKNESS=0x28 # short
    OFFSET_MODULUS = 0x2A # float
    OFFSET_DIAMETER = 0x2E # short
    OFFSET_ANGLE = 0x30 # short
    OFFSET_X_START = 0x32 # float
    OFFSET_X_STOP = 0x36 # float
    OFFSET_X_STEP = 0x3A # float
    OFFSET_TIME = 0x3E # float
    OFFSET_FIRST = 0x42 # byte
    OFFSET_WORKINGDISTANCE = 0x43 # float
    OFFSET_BOW = 0x47 # float
    OFFSET_MYSTERY_F2 = 0x4B # float
    OFFSET_LOWINTENSITYALARM = 0x4F # float
    OFFSET_AMINUSB = 0x53 # 50 floats
    OFFSET_APLUSB = OFFSET_AMINUSB + 50*4 # 50 floats
    OFFSET_MYSTERY_F3 = 0x1E3
    OFFSET_MYSTERY_F4 = 0x1E7
    MAGIC_SCALER = 0.04517943889
    MAGIC_OFFSET = 0.0001789666667

    TimeOffset = 3*3600 # There seems to be a 3 hour time offset

    def __init__(self, sampleId = None, comment = None):
        self.sampleId = sampleId
        self.comment = comment
        self.maxBow = None
        self.radius = 0
        self.temperature = None
        self.filmThickness = 0
        self.waferThickness = 0
        self.modulus = None
        self.angle = 0
        self.diameter = 0
        self.xStart = 0
        self.xStop = 0
        self.xStep = 0
        self.y = []
        self.intensity = []
        self._bow = None

    @property
    def x(self):
        '''Return the scan points [m]'''
        #xr = np.arange(self.xStart, self.xStop, self.xStep)
        xr = np.linspace(self.xStart, self.xStop, len(self.y))
        #print len(xr), len(self.y)
        assert(len(xr)==len(self.y))
        return xr

    @staticmethod
    def _unpackFloats(data, position, n=1):
        if n==1:
            return struct.unpack('<f', data[position:position+4])[0]
        else:
            return np.array(struct.unpack('<%df' % n, data[position:position+n*4]))

    @staticmethod
    def _unpackShorts(data, position, n=1):
        if n==1:
            return struct.unpack('<H', data[position:position+2])[0]
        else:
            return np.array(struct.unpack('<%dH' % n, data[position:position+n*2]))

    @staticmethod
    def _unpackInts(data, position, n=1):
        if n==1:
            return struct.unpack('<I', data[position:position+4])[0]
        else:
            return np.array(struct.unpack('<%dI' % n, data[position:position+n*4]))

    def decodeFromBinary(self, data):
        '''Decode the scan record from a binary blob (sca file)'''
        assert(len(data) == Scan.RECORD_LENGTH)
        i = min(data.index('\0',0), Scan.OFFSET_COMMENT)
        self.sampleId = data[0:i]
        i = data.index('\0', Scan.OFFSET_COMMENT)
        self.comment = data[Scan.OFFSET_COMMENT:i]
        self.mysteryS = Scan._unpackShorts(data, Scan.OFFSET_MYSTERY_S)
        self.first = bool(struct.unpack('<B', data[Scan.OFFSET_FIRST])[0])
        self.maxBow = Scan._unpackFloats(data, Scan.OFFSET_BOW)*1E-6
        self.time = time.localtime(Scan._unpackInts(data, Scan.OFFSET_TIME)-self.TimeOffset)
        self.mysteryF2 = Scan._unpackFloats(data, Scan.OFFSET_MYSTERY_F2)
        self.mysteryF3 = Scan._unpackFloats(data, Scan.OFFSET_MYSTERY_F3)
        self.mysteryF4 = Scan._unpackFloats(data, Scan.OFFSET_MYSTERY_F4)
        self.radius = Scan._unpackFloats(data, Scan.OFFSET_RADIUS)
        self.lowIntensityAlarm = Scan._unpackFloats(data, Scan.OFFSET_LOWINTENSITYALARM)
        self.temperature = Scan._unpackShorts(data, Scan.OFFSET_TEMPERATURE) + 273.15 # From C to Kelvin
#        self.filmThickness = Scan._unpackShorts(data, Scan.OFFSET_FILMTHICKNESS)*1E-10 # From Angstrom to meter
        self.filmThickness = Scan._unpackInts(data, Scan.OFFSET_FILMTHICKNESS)*1E-10 # From Angstrom to meter
        self.workingDistance = Scan._unpackFloats(data, Scan.OFFSET_WORKINGDISTANCE) * 0.0254 # From inch to meter
        self.waferThickness = Scan._unpackShorts(data, Scan.OFFSET_WAFERTHICKNESS) * 1E-6 # From micron to meter
        self.detectorDistance = self.workingDistance+self.waferThickness
        self.modulus = Scan._unpackFloats(data, Scan.OFFSET_MODULUS)*1E11 # To MPa
        self.angle = Scan._unpackShorts(data, Scan.OFFSET_ANGLE)
        self.diameter = Scan._unpackShorts(data, Scan.OFFSET_DIAMETER) * 1E-3
        self.xStart = Scan._unpackFloats(data, Scan.OFFSET_X_START) * 1E-3
        self.xStop = Scan._unpackFloats(data, Scan.OFFSET_X_STOP) * 1E-3
        self.xStep = Scan._unpackFloats(data, Scan.OFFSET_X_STEP) * 1E-2 # Not sure about this one
        self.AmB = Scan._unpackFloats(data, Scan.OFFSET_AMINUSB, 50)
        self.ApB = Scan._unpackFloats(data, Scan.OFFSET_APLUSB, 50)
        self.intensity = self.ApB
        self.shift = self.AmB / self.ApB # (A-B) / (A+B)
        self.y = self.shift
        self.calibrateDeflection()
        #self.deflection = self.shift*self.MAGIC_SCALER + self.MAGIC_OFFSET
        self.deflection = self.scaleFactor * self.shift + self.MAGIC_OFFSET

    def calibrateDeflection(self):
        i = self.intensity > self.lowIntensityAlarm
        assert(np.count_nonzero(i) > 3)
        fit = np.polyfit(self.x[i], self.shift[i], 1)
        slope = fit[0]
        R = 1/(slope)
        self.scaleFactor = R/self.radius
        print "Scale factor:" , self.scaleFactor

    @property
    def bow(self):
        '''Calculates bow as a function of x from deflection data'''
        if self._bow is None:
            self._bow = integrate.cumtrapz(self.deflection-np.mean(self.deflection), x=self.x)
            self._bow = np.insert(self._bow, 0, 0)
        return self._bow


    def table(self):
        '''Return data in a table format'''
        s  = '%16s\t' % self.sampleId
        s += '%12s\t' % self.comment
        s += time.strftime('%Y-%m-%d %H:%M:%S\t', self.time)
        s += "%d\t" % (self.temperature-273.15)
        s += "%d\t" % self.angle
        s += "%d\t" % (self.waferThickness*1E6)
        s += "%e\t" % self.modulus
        s += "%d\t" % (self.filmThickness*1E10)
        s += "%g\t" % (self.diameter*1E3)
        s += "%g\t" % (self.xStart*1E3)
        s += "%g\t" % (self.xStop*1E3)
        s += "%.3f\t" % (self.maxBow*1E6)
        s += "%.3f\t" % self.radius
        s += "%f\t" % self.lowIntensityAlarm
        s += "%f\t" % self.workingDistance
        s += "%f\t" % self.mysteryF2
        s += "%f\t" % self.mysteryF3
        s += "%f" % self.mysteryF4
        return s

    @staticmethod
    def tableHeader():
        '''Return a header for the table.'''
        s  = 'Sample ID\t'
        s += 'Comment\t'
        s += 'Time\t'
        s += 'Temperature [C]\t'
        s += 'Angle [deg]\t'
        s += 't_Wafer [um]\t'
        s += 'Modulus [Pa]\t'
        s += 't_Film [A]\t'
        s += 'Diameter [mm]\t'
        s += 'Scan start [mm]\t'
        s += 'Scan stop [mm]\t'
        s += 'Bow [um]\t'
        s += 'Radius [m]\t'
        s += 'Low intensity alarm [V]\t'
        s += 'Working distance [m]\t'
        s += 'Mystery F2\t'
        s += 'Mystery F3\t'
        s += 'Mystery F4\t'
        return s

    def __str__(self):
        s  = "Sample ID: %s" % self.sampleId
        s += "\tComment: %s" % self.comment
        s += "\tTime: %s" % time.strftime('%Y-%m-%d %H:%M:%S', self.time)
        s += "\tTemperature: %d C" % (self.temperature-273.15)
        s += "\tAngle: %d deg" % self.angle
        s += "\tt_Wafer: %d um" % (self.waferThickness*1E6)
        s += "\tModulus: %e Pa" % self.modulus
        s += "\tt_Film: %d A" % (self.filmThickness*1E10)
        s += "\tDiameter:%g mm" % (self.diameter*1E3)
        s += "\tScan range: %g to %g (step %g) mm" % (self.xStart*1E3,self.xStop*1E3,self.xStep*1E3)
        s += "\tBow: %.3f um" % (self.bow*1E6)
        s += "\tRadius: %.3f m" % self.radius
        s += "\tLow Intensity alarm: %f V" % self.lowIntensityAlarm
        s += "\tWorking distance: %f m" % self.workingDistance
        s += "\tFirst: %s" % str(self.first)
        s += "\tMysteryS: %f" % self.mysteryS
        s += "\tMystery F2: %f" % self.mysteryF2
        s += "\tMystery F3: %f" % self.mysteryF3
        s += "\tMystery F4: %f" % self.mysteryF4
        return s

    def __repr__(self):
        return str(self)

class ScanImporter(object):
    '''Import scan records from Tencor .sca file.'''
    OFFSET=0x64

    def __init__(self, fileName):
        self.scans = []
        self.sampleIds = []
        self.readFile(fileName)

    def appendScan(self, scan):
        self.scans.append(scan)
        self.sampleIds.append(scan.sampleId)

    def readFile(self, fileName):
        with open(fileName, 'r') as f:
            d = f.read()
        l = len(d)
        nRecords = (l-ScanImporter.OFFSET) / Scan.RECORD_LENGTH
        lExpected = nRecords*Scan.RECORD_LENGTH + ScanImporter.OFFSET
        magic = struct.unpack('<I', d[0:4])[0]
        if magic != 0x00C8A5A5:
            print "Warning: Incorrect magic bytes in header!"
        if l < lExpected:
            print "Warning: File seems to be too short. Missing %d bytes." % lExpected-l
        elif l > lExpected:
            print "Warning: File seems to be too long. Have %d extra bytes." % l-lExpected
        for n in range(nRecords):
            data = d[ScanImporter.OFFSET+Scan.RECORD_LENGTH*n:ScanImporter.OFFSET+Scan.RECORD_LENGTH*(n+1)]
            scan = Scan()
            scan.decodeFromBinary(data)
            self.appendScan(scan)

    def numberOfScans(self):
        return len(self.scans)

    def __len__(self):
        return self.numberOfScans()

    def __getitem__(self, item):
        if type(item) is int:
            return self.scans[item]
        elif type(item) is str:
            r = []
            for i,sampleId in enumerate(self.sampleIds):
                if sampleId == item:
                  r.append(self.scans[i])
            return r
        else:
            raise KeyError('Unsupported scan identifier: %s' % str(item))

def importTxt(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    x = []
    y = []
    for l in lines[2:]:
        i = l.find(':')
        r = l[i+1:].split()
        x.append(float(r[0]))
        y.append(float(r[1]))
    return np.array(x), np.array(y)

if __name__ == '__main__':
    import matplotlib.pyplot as mpl
    wafer = '2686_J25'
    fileName = 'KLK/%s.DAT' % wafer
    t = TencorTable(fileName)
    t.saveToExcel('%s.xls' % wafer)


#    case = 2
#    if case == 1:
#        fn = 'FJ/FELIX.SCN'
#        scans = ScanImporter(fn)
#        sampleIds = ['2314 1090', '2314 2080', '2314 3070']
#        scans = [scan for scan in scans if scan.sampleId in sampleIds]
#        textFiles = ['FJ/%s.TXT' % (''.join(sampleId.split()) ) for sampleId in sampleIds]
#    elif case == 2:
#        fn = 'FJ2/FJ2.SCN'
#        scans = ScanImporter(fn)
#        print "Number of scans in file:", scans.numberOfScans()
#        textFiles = ['FJ2/%s.TXT' % scan.sampleId for scan in scans]
#    print Scan.tableHeader()
#    for i, scan in enumerate(scans):
#        print scan.table()
#    for i,scan in enumerate(scans):
#
#        #run = scan.sampleId
#        #fn = 'FJ2/%s.TXT' % run
#        fn = textFiles[i]
#        #print "Text file:", fn
#        xText, yText = importTxt(fn)
#        f = np.polyfit(scan.y, yText, 1)
#        #alpha = scan.mysteryF4/10*scan.workingDistance*scan.y - scan.mysteryF3/10
#        #f = linearFit(scan.y, yText)
#        #f = [f,0]
#        print "%s\t%s\t%.7f\t%.7f" % (scan.sampleId, scan.comment, f[0], f[1])
#        mpl.subplot(311)
#        mpl.title('Text file %s, scan %s (%s)' % (fn, scan.sampleId, scan.comment))
#        mpl.plot(scan.x*1E3, scan.y, '.-')
#        mpl.ylabel('y$_{Scan}$ [?]')
#        mpl.subplot(312)
#        mpl.plot(xText, yText, '.-', label='text')
#        mpl.plot(scan.x*1E3, f[0]*scan.y+f[1],'.-', label='re-scaled scan')
#        mpl.plot(scan.x*1E3, scan.deflection,'.-', label='alpha')
#        mpl.ylabel('y$_{Text}$ [mrad]')
#        mpl.legend(loc='lower right')
#        mpl.xlabel('x [mm]')
#        mpl.subplot(313)
#        mpl.plot(scan.y, yText, 'o')
#        mpl.plot(scan.y, np.polyval(f, scan.y), '-', label='fit=%f*y_{Scan} + %f' % (f[0],f[1]))
#        mpl.xlabel('y$_{Scan}$ [?]')
#        mpl.ylabel('y$_{Text}$ [mrad]')
#        mpl.legend(loc='lower right')
#        mpl.show()
#        mpl.close()
