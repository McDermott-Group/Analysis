# Import Plot.ly functions for both online/offline
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

# Sign in to Plotly for API access (required if publishing online)

# Import extra modules as needed
import numpy as np
import matplotlib as mpl
import scipy as sp

import re
import os
import pandas

from WaferMap import WaferMap
from dataChest import dataChest

init_notebook_mode()

class probeTest(object):
    
    def __init__(self, path, bounds=None):
        """bounds should be in the form {area: (lower,upper)}"""
        chest = dataChest(path[:-1])
        chest.openDataset(path[-1])
        # data should be in format (die, area, range, resistance)
        self.data = [(row[3],row[1],row[0],str(row[2])) for row in chest.getData()
                     if not np.isnan(row[2]) \
                     and (bounds==None \
                          or row[1] not in bounds \
                          or (row[3]>bounds[row[1]][0] \
                              and row[3]<bounds[row[1]][1]))]
        # data in format (resistance, area, die, range)
        self.odd = True
        self.inner_diameter = 65
        self.pitchX = 6.2
        self.pitchY = 6.2
        try:
            self.odd = chest.getParameter('Odd')
            self.inner_diameter = chest.getParameter('Inner Diameter')
            self.pitchX = chest.getParameter('Pitch X')
            self.pitchY = chest.getParameter('Pitch Y')
        except:
            pass
        # what areas do we have?
        self.areas = []
        for row in self.data:
            a = float(row[1])
            if a not in self.areas:
                self.areas.append(a)
        
    def calcMean(self):
        """Returns a dictionary of {area:average}"""
        avg = {a:{'sum':0, 'n':0} for a in self.areas}
        for row in self.data:
            r,a = float(row[0]), float(row[1])
            avg[a]['sum'] += r
            avg[a]['n'] += 1
        for a in avg:
            avg[a] = avg[a]['sum']/avg[a]['n']
        return avg
    
    def theoreticalR(self, Jc, gap=380):
        """Calculate what we expect theoretically for room temperature R,
        based on area (in um^2), gap 2\Delta/e (in uV), and critical current
        density (in A/cm^2)."""
        return {a:np.pi/4*gap*1e-6/(Jc*1e4*a*1e-12) for a in self.areas}
    
    def resistancePlot(self):
        x = [1/float(line[1]) for line in self.data]
        y = [float(line[0]) for line in self.data]
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
                mode = 'line',
                name = 'Fit'
            )
        fig = dict(data = [scatter,fit], 
                   layout = {'title':'',
                            'xaxis':{'title':'1/JJ Area [1/um^2]'},
                            'yaxis':{'title':'Room Temperature Resistance [\Omega]'}
                            }
                  )
        iplot(fig)
    
    def resistanceHistogram(self, area=None, binWidth=10):
        if area:
            areaList = [area]
        else:
            areaList = self.areas
        hist = [go.Histogram(
            x = [float(line[0]) for line in self.data if line[1]==a],
            autobinx=False,
            xbins=dict(
                start=0,
                end=1000,
                size=binWidth
            )
        ) for a in areaList]
        iplot(hist)
        
    def dieMap(self, area, fn=np.mean):
        w = WaferMap(inner_diameter=self.inner_diameter, pitch=self.pitchX, odd=self.odd, title=None)
        data = [row for row in self.data if float(row[1])==float(area)]
        for row in data:
            die = row[2]
            dieData = [row[0] for row in data if row[2] == die]
            w[die] = fn(dieData)# - self.calcMean()[area]
        w.show()
