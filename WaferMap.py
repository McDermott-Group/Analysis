import numpy as np

# Import Plot.ly functions for both online/offline
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

class WaferMap(object):
    def __init__(self, wafer_diameter=76.2, sputter_diameter=71.12, inner_diameter=65., pitch=6.2, odd=True, title=None):
        """Old Label is where the second index in C3 means 3 down from the top die, not 3 down from the top of the grid."""
        self.wafer_diameter = wafer_diameter
        self.sputter_diameter = sputter_diameter
        self.inner_diameter = inner_diameter
        self.pitch = pitch
        self.title = title
        if odd:
            self.ndie = int(np.floor(inner_diameter/pitch))
            self.center_die = np.floor(self.ndie/2.)
        else:
            self.ndie = int(2.*np.floor((inner_diameter-pitch)/2/pitch))
            self.center_die = self.ndie/2.-0.5
        self.dieMapMask = np.ones((self.ndie+2,self.ndie+2), dtype=bool) # is die on wafer?
        self.dieMapData = np.zeros((self.ndie+2,self.ndie+2)) # actual data stored here
        for x in range(self.ndie+2):
            for y in range(self.ndie+2):
                x_cent = (x-1-self.center_die)*pitch
                y_cent = (y-1-self.center_die)*pitch
                if np.sqrt(x_cent**2 + y_cent**2) >= inner_diameter/2.:
                    self.dieMapMask[y][x] = False
                    self.dieMapData[y][x] = np.nan
                else:
                    self.dieMapData[y][x] = np.nan
    def __setitem__(self, key, value):
        colIndex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.find(key.upper()[0])
        rowIndex = int(key[1:]) - 1
        self.dieMapData[rowIndex+1][colIndex] = value
        return
        if self.dieMapMask[colIndex][rowIndex+1]:
            self.dieMapData[rowIndex+1][colIndex] = value
        else: 
            raise KeyError('die not on wafer: '+key)
    def __getitem__(self, key):
        colIndex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.find(key.upper()[0])
        rowIndex = int(key[1:]) - 1
        if self.dieMapMask[colIndex][rowIndex+1]:
            return self.dieMapData[rowIndex+1][colIndex]
        else: 
            raise KeyError('die not on wafer')
    def show(self):
        annotations = []
        pitch = self.pitch
        center_die = self.center_die
        wafer_diameter = self.wafer_diameter
        sputter_diameter = self.sputter_diameter
        trace = go.Heatmap(
                z = self.dieMapData,
                x = [c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'],
                # y = [str(self.ndie-i)+'a' for i in range(self.ndie)]
            )
        if self.title is None: 
            self.title = 'Wafer Die Map'
        fig = dict(data = [trace], 
                   layout = {'title': self.title,
                             'width':500,
                             'height':500,
                             'autosize':False,
                             'xaxis': {'ticks':'', 'side':'top'},
                             'yaxis':{'autorange':'reversed'},
                             'annotations':annotations,
                             'shapes': [
                                # outer circle
                                {
                                    'type': 'circle',
                                    'xref': 'x',
                                    'yref': 'y',
                                    'x0': center_die-wafer_diameter/2./pitch+1,
                                    'y0': center_die-wafer_diameter/2./pitch+1,
                                    'x1': center_die+wafer_diameter/2./pitch+1,
                                    'y1': center_die+wafer_diameter/2./pitch+1,
                                    'line': {
                                        'color': 'rgba(50, 171, 96, 1)',
                                    },
                                },
                                # middle circle
                                {
                                    'type': 'circle',
                                    'xref': 'x',
                                    'yref': 'y',
                                    'x0': center_die-sputter_diameter/2./pitch+1,
                                    'y0': center_die-sputter_diameter/2./pitch+1,
                                    'x1': center_die+sputter_diameter/2./pitch+1,
                                    'y1': center_die+sputter_diameter/2./pitch+1,
                                    'line': {
                                        'color': 'rgba(171, 50, 96, 1)',
                                    },
                                },]
                            }
                  )
        iplot(fig)