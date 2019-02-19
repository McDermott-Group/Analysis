# Import Plot.ly functions for both online/offline
import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.graph_objs as go
import plotly.tools as tls
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

#import matplotlib.pyplot as plt
#%matplotlib inline

# Sign in to Plotly for API access (required if publishing online)
# Comment out to prevent accidental publishing
# py.sign_in('iamed18','o9sd2728op')

# Import extra modules as needed
import numpy as np
import matplotlib as mpl
import scipy as sp

# import peakutils
import re
import os
import pandas

from WaferMap import WaferMap
from dataChest import dataChest

init_notebook_mode()




class probeTest(object):
    
    def __init__(self, path, bounds={}):
        """bounds should be in the form {area: (lower,upper)}"""
        chest = dataChest(path[:-1])
        chest.openDataset(path[-1])
        dataOrder = ('die', 'index', 'area', 'range', 'resistance')
        self.data = [dict(zip(dataOrder,row)) for row in chest.getData()
                     if not np.isnan(row[2]) \
                     and (row[2] not in bounds \
                          or (row[4]>bounds[row[2]][0] \
                              and row[4]<bounds[row[2]][1]))]
        # data in format (resistance, area, die, range)
        self.odd = False
        self.inner_diameter = 65
        self.pitchX = 6.25
        self.pitchY = 6.25
        try:
            self.odd = chest.getParameter('Odd')
            self.inner_diameter = chest.getParameter('Inner Diameter')
            self.pitchX = chest.getParameter('Pitch X')
            self.pitchY = chest.getParameter('Pitch Y')
        except:
            pass
    
    def removeOverwrittenData(self):
        for i in range(len(self.data)-1,-1,-1):
            for j in range(i-1,-1,-1):
                if (self.data[i]['die'] == self.data[j]['die']
                        and self.data[i]['index'] == self.data[j]['index']):
                    del self.data[j]
                    break
    
    def addSubstrateCalibration(self, path):
        """This should be a path to a normal probe file that only has a single die with a measurement at each index."""
        chest = dataChest(path[:-1])
        chest.openDataset(path[-1])
        dataOrder = ('die', 'index', 'area', 'range', 'resistance')
        self.calibration_data = [dict(zip(dataOrder,row)) for row in chest.getData()]
        cal_dict = {row['index']:row['resistance'] for row in self.calibration_data}
        for row in self.data:
            row['resistance'] = 1./( 1./row['resistance'] + 1./cal_dict[row['index']] )
        
    def calcMean(self):
        """Returns a dictionary of {area:average}"""
        areas = []
        for row in self.data:
            if row['area'] not in areas:
                areas.append(row['area'])
        avg = {a:{'sum':0, 'n':0} for a in areas}
        for row in self.data:
            a = float(row['area'])
            avg[a]['sum'] += float(row['resistance'])
            avg[a]['n'] += 1
        for a in avg:
            avg[a] = avg[a]['sum']/avg[a]['n']
        return avg
    
    def theoreticalR(self, Jc, gap=380):
        """Calculate what we expect theoretically for room temperature R,
        based on area (in um^2), gap 2\Delta/e (in uV), and critical current
        density (in A/cm^2)."""
        areas = []
        for row in self.data:
            if row['area'] not in areas:
                areas.append(row['area'])
        return {a:np.pi/4*gap*1e-6/(Jc*1e4*a*1e-12) for a in areas}
    
    def resistancePlot(self):
        """Old version"""
        x = [1/float(line['area']) for line in self.data]
        y = [float(line['resistance']) for line in self.data]
        scatter = go.Scatter(
                x = x,
                y = y,
                mode = 'markers',
                name = 'Data'
            )
        xfit = np.unique(x)
        m, b = np.polyfit(x, y, 1)
        yfit = [m*xi+b for xi in xfit]
        fit = go.Scatter(
                x = xfit,
                y = yfit,
                mode = 'lines',
                name = 'Fit'
            )
        fig = dict(data = [scatter,fit], 
                   layout = {'title':'',
                            'xaxis':{'title':'1/JJ Area [1/um^2]'},
                            'yaxis':{'title':'Room Temperature Resistance [\Omega]'}
                            }
                  )
        iplot(fig)
    
    def resistancePlot(self):
        """This just allows grouping by index (or die, etc)"""
        scatter = [go.Scatter(
                x = [float(line['area']) for line in self.data if line['index']==die],
                y = [1/float(line['resistance']) for line in self.data if line['index']==die],
                mode = 'markers',
                name = 'Data'
            ) for die in np.unique([row['index'] for row in self.data])]
        
        # linear fit based on all data
        x = [float(line['area']) for line in self.data]
        y = [1/float(line['resistance']) for line in self.data]
        xfit = np.unique(x)
        m, b = np.polyfit(x, y, 1)
        yfit = [m*xi+b for xi in xfit]
        fit = go.Scatter(
                x = xfit,
                y = yfit,
                mode = 'lines',
                name = 'Fit'
            )
        fig = dict(data = scatter+[fit], 
                   layout = {'title':'',
                            'xaxis':{'title':'Junction Area [um^2]'},
                            'yaxis':{'title':'1/Room Temperature Resistance [$1/\Omega$]'}
                            }
                  )
        iplot(fig)
    
    def displayDataTable(self):
        colors = ['rgba(153,{},255,100)'.format( 153+50*(ord(row['die'][0])%2) ) for row in self.data]
        table = go.Table(
            header=dict(values=['Die', 'Area', 'Resistance']),
            cells=dict(values=[[row['die'] for row in self.data],
                               [row['area'] for row in self.data],
                               ['{:.2f}k'.format(row['resistance']/1000) for row in self.data]],
                       fill = dict(color = [colors, colors, colors])
                      ))
        iplot([table], filename = 'basic_table')
    
    def resistanceHistogram(self, area=None, binWidth=1000):
        if area:
            areaList = [area]
        else:
            areaList = []
            for row in self.data:
                if row['area'] not in areaList:
                    areaList.append(row['area'])
        hist = [go.Histogram(
            x = [float(line['resistance']) for line in self.data if line['area']==a],
            autobinx=False,
            xbins=dict(
            #    start=0,
            #    end=1000,
                size=binWidth
            )
        ) for a in areaList]
        iplot(hist)
        
    def dieMap(self, area, fn=np.mean):
        w = WaferMap(inner_diameter=self.inner_diameter, pitch=self.pitchX, odd=self.odd, title=None)
        dies = np.unique([row['die'] for row in self.data])
        for die in dies:
            resList = [row['resistance'] for row in self.data if row['die'] == die and float(row['area'])==float(area)]
            w[die] = fn(resList)# - self.calcMean()[area]
        w.show()