import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
import matplotlib.pyplot as plt
import noiselib
from numba import jit
import time


class ImpactEvent(object):
    
    def __init__(self, pdfs, track, energy, qubit_xys, charge_map):
    
        self.PDFs = pdfs
        self.track = track
        self.energy = energy
        self.qubit_xys = qubit_xys
        self.charge_map = charge_map
    
    def set_track(self, track, energy):
    
        self.track = track
        self.energy = energy
        if hasattr(self, 'charge'):
            del self.charge
    
    def get_track(self):
    
        c = self.energy
        colorFloats = np.interp(c, (c.min(), c.max()), (0,1))
        color = cm.autumn(colorFloats)
        return self.track, color
    
    def get_charge(self, fQ=1., cache=False):
        """ Generates and returns positions and polarities of charges from track 
        and energy. """
    
        if hasattr(self, 'charge'):
            return self.charge
        xyz = []
        polarity = []
        for i in range(len(self.track)):
            x,y,z = self.track[i]
            res = self.diffuseCharge(x,y,z,self.energy[i],fQ=fQ,epseh=3.6,verbose=False)
            polarity += [np.full(res['electrons'].shape[1], -1.), 
                         np.full(res['holes'].shape[1], 1.) ]
            xyz += [res['electrons'].T, res['holes'].T]
        xyz = np.concatenate(xyz)
        polarity = np.concatenate(polarity)
        if cache:
            self.charge = xyz, polarity
            print( xyz.nbytes/1e6 )
        return xyz, polarity
    
    def get_charge_on_qubits(self, charge_xyz, charge_polarity):
        """ Return the charge that lands directly on each qubit island. """
    
        q = []
        for i,xy_q in enumerate(self.qubit_xys):
            dx,dy,dz = (charge_xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
            r = np.hypot(dx,dy)
            qs = (r < 0.07) & (dz == 0)
            q.append( (np.sum(qs[charge_polarity<0]),
                       np.sum(qs[charge_polarity>0])) )
        return q
    
    def get_induced_charge_on_qubits(self, charge_xyz, charge_polarity, sum=True):
        """ Calculate the induced offset charge on each qubit from the distribution
        charge_xyz and polarity charge_polarity of the charges.  If sum==True, 
        the induced charge will be summed before being returned, giving the total
        offset charge.  If False, an array of induced charges will be returned
        for each qubit. """
    
        q = []
        for i,xy_q in enumerate(self.qubit_xys):
            dx,dy,dz = (charge_xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
            r = np.hypot(dx,dy)
            qs = self.charge_map(r,-dz)
            if sum:
                q.append( np.sum(qs*charge_polarity) )
            else:
                q.append(qs)
        return np.array(q)
    
    def getPDF(self, z, charge_type='electrons'):
        
        zcuts = self.PDFs['zarray']
        zind = np.argmin(np.abs(z-zcuts))
        return self.PDFs[zcuts[zind]][charge_type]
    
    def diffuseCharge(self, x, y, z, energy, fQ = 1.0, epseh=3.6, verbose=True):
    
        PDFs = self.PDFs
        
        neh = int(np.floor(energy/epseh)*fQ)
        if(verbose):
            print('Event Energy: '+str(energy)+' eV')
            print('Electron-Hole Pairs: '+str(neh))
        
        results=dict()
        for charge_type in ['electrons','holes']:
            Hn=self.getPDF(z,charge_type=charge_type)
            Hr = Hn.ravel()
            # vals = np.random.choice(Hr.size, p=Hr.astype('float64'), size=neh)
            r = np.random.random(neh)
            cumsum = np.cumsum(Hr)
            vals = np.searchsorted(cumsum/np.sum(Hr), r)
            inds = np.unravel_index(vals,Hn.shape)

            sim_distx = PDFs['x'][inds[0]]
            sim_disty = PDFs['y'][inds[1]]
            sim_distz = PDFs['z'][inds[2]]
        
            xf = x+sim_distx
            yf = y+sim_disty
            zf = sim_distz
            
            xf += +np.random.rand(len(xf))*PDFs['dx'] - PDFs['dx']/2.0
            yf += +np.random.rand(len(yf))*PDFs['dy'] - PDFs['dy']/2.0
            zf += +np.random.rand(len(zf))*PDFs['dz'] - PDFs['dz']/2.0
            
            zmax = self.PDFs['thickness']/2.
            zmin = -self.PDFs['thickness']/2.
            xf[zf > zmax] = (x+(zf-z)*(xf-x))[zf > zmax]
            xf[zf < zmin] = (x+(zf-z)*(xf-x))[zf < zmin]
            yf[zf > zmax] = (y+(zf-z)*(yf-y))[zf > zmax]
            yf[zf < zmin] = (y+(zf-z)*(yf-y))[zf < zmin]
            zf[zf > zmax] = zmax
            zf[zf < zmin] = zmin
            
            if charge_type == 'electrons':
                results['electrons']=np.array([xf,yf,zf])
            else:
                results['holes']=np.array([xf,yf,zf])
        
        return results


class MainWindow(gl.GLViewWidget):
    
    def __init__(self, model, qubit_xys):
    
        super(MainWindow, self).__init__()
        # self.setGeometry(300, 300, 250, 150)
        self.show()
        self.model = model
        
        # g = gl.GLGridItem()
        # g.setSize(10,10,10)
        # self.addItem(g)
        
        self.track = gl.GLScatterPlotItem( size=5 )
        self.track.setGLOptions('translucent')
        self.addItem(self.track)
        
        self.charge = gl.GLScatterPlotItem( size=2 )
        self.charge.setGLOptions('translucent')
        self.addItem(self.charge)
        
        ax = gl.GLAxisItem()
        ax.setSize(10,10,10)
        self.addItem(ax)
        
        self.draw_chip_lines()
        
        for x,y in qubit_xys:
            self.draw_qubit(x,y)
    
    def draw_chip_faces(self, h=0.375/2., L=6.25/2.):
        
        vertexes = np.array([[-L,-L,-h],
                             [-L,-L, h],
                             [-L, L,-h],
                             [-L, L, h],
                             [ L,-L,-h],
                             [ L,-L, h],
                             [ L, L,-h],
                             [ L, L, h]])
        faces = np.array([[0,2,4], [2,4,6],
                          [0,1,5], [0,4,5],
                          [0,2,3], [0,1,3]])
        return gl.GLMeshItem(vertexes=vertexes, faces=faces)
    
    def draw_chip_lines(self, h=0.375/2., L=6.25/2.):
        vertexes = np.array([[-L,-L,-h],
                             [ L,-L,-h],
                             [ L, L,-h],
                             [-L, L,-h],
                             [-L,-L,-h],
                             [np.nan,np.nan,np.nan],
                             [-L,-L, h],
                             [ L,-L, h],
                             [ L, L, h],
                             [-L, L, h],
                             [-L,-L, h],
                             [np.nan,np.nan,np.nan],
                             [-L,-L,-h],
                             [-L,-L, h],
                             [np.nan,np.nan,np.nan],
                             [-L, L,-h],
                             [-L, L, h],
                             [np.nan,np.nan,np.nan],
                             [ L,-L,-h],
                             [ L,-L, h],
                             [np.nan,np.nan,np.nan],
                             [ L, L,-h],
                             [ L, L, h]])
        self.addItem( gl.GLLinePlotItem(pos=vertexes) )
    
    def draw_qubit(self, x, y, r_inner=0.07, r_outer=0.0905, h=0.375/2.):
        
        outer_circle = self._circle((x,y), r_outer)
        outer_circle = np.concatenate([outer_circle, [outer_circle[0]]])
        self.addItem( gl.GLLinePlotItem(pos=outer_circle) )
        
        inner_circle = self._circle((x,y), r_inner)
        vertexes = np.concatenate([[[x,y,h]], inner_circle])
        pts = range(1,len(vertexes))
        faces = np.array([ [0,pts[i],pts[i+1]] for i in range(-1,len(pts)-1) ])
        self.addItem( gl.GLMeshItem(vertexes=vertexes, faces=faces) )
    
    def _circle(self, xy, r, h=0.375/2., n=100):
        """Returns points for a circle."""
        x = xy[0] + r*np.cos(np.linspace(0,2*np.pi,n))
        y = xy[1] + r*np.sin(np.linspace(0,2*np.pi,n))
        h = np.full(n, h)
        return np.transpose([x,y,h])
        
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Left:
            self.model.show_previous()
        elif event.key() == QtCore.Qt.Key_Right:
            self.model.show_next()
        else:
            return
    
    def plot_track(self, track_xyz, track_color, charge_xyz, charge_polarity):
    
        self.track.setData(pos=track_xyz,
                           color=track_color)
        self.charge.setData(pos=charge_xyz,
                            color=cm.PiYG(charge_polarity))


class Controller(QtGui.QApplication):
    
    def __init__(self, sys_argv, plot=False, calc_all=True,
                 pdfs_file='sim_data/ChargePDFs_150.npy',
                 event_files=['sim_data/Gamma.txt', 'sim_data/Gamma_10deg.txt'],
                 fQ=1.):
                 
        super(Controller, self).__init__(sys_argv)
    
        self.qubit_xys = [( 1.564,  0.570),
                          ( 1.564, -0.070),
                          (-1.564, -0.080),
                          (-1.564, -0.420)]
        
        self.fQ = fQ
        
        PDFs = np.load(pdfs_file, allow_pickle=True)
        PDFs = PDFs.tolist()
        
        q_map = noiselib.loadmat('sim_data/charge_map.mat')
        q = q_map['charge_mat']
        q[~np.isfinite(q)] = -1. # right below center of qubit
        def map_charge(r, z, q, rp, zp): 
            # my own bilinear interpolation, much faster
            # https://math.stackexchange.com/questions/3230376/
            #   interpolate-between-4-points-on-a-2d-plane
            ri = np.searchsorted(r, rp)
            zi = np.searchsorted(z, zp)
            dr, dz = r[1]-r[0], z[1]-z[0]
            rpp = (rp-r[np.clip(ri-1,0,r.size-1)])/dr
            zpp = (zp-z[np.clip(zi-1,0,z.size-1)])/dz
            q1 = q[np.clip(zi-1,0,z.size-1), np.clip(ri-1,0,r.size-1)]
            q2 = q[np.clip(zi-1,0,z.size-1), np.clip(ri  ,0,r.size-1)]
            q3 = q[np.clip(zi  ,0,z.size-1), np.clip(ri  ,0,r.size-1)]
            q4 = q[np.clip(zi  ,0,z.size-1), np.clip(ri-1,0,r.size-1)]
            return (1-rpp)*(1-zpp)*q1 + rpp*(1-zpp)*q2 + (1-rpp)*zpp*q3 + rpp*zpp*q4
        charge_map = lambda rp,zp: map_charge( q_map['r']/1000., 
                                               q_map['z']/1000., 
                                               q, rp, zp )
        
        data = []
        for f in event_files:
            data.append( np.loadtxt(f, skiprows=1) )
        data = np.concatenate(data)
        self.split_data = np.split(data, np.where(np.diff(data[:,0]))[0]+1)
        self.events = []
        for event in self.split_data:
            e = ImpactEvent(PDFs, event[:,[2,3,4]], event[:,5] * 1e6,
                            self.qubit_xys, charge_map)
            self.events.append(e)
        
        self.i = -1
        
        if plot:
            self.view = MainWindow(self, self.qubit_xys)
            self.show_next()
        
        if calc_all:
            self.e = ImpactEvent(PDFs, 1, 1, self.qubit_xys, charge_map)
            self.thread = threading.Thread(target=self.analyze_all_events)
            self.thread.start()
    
    def show_previous(self):
        
        self.i -= 1
        self.show_event()
    
    def show_next(self):
        
        self.i += 1
        self.show_event()
        
    def show_event(self):
    
        event = self.events[ self.i % len(self.split_data) ]
        track_xyz, track_color = event.get_track()
        charge_xyz, charge_polarity = event.get_charge(fQ=self.fQ)
        self.view.plot_track(track_xyz, track_color, charge_xyz, charge_polarity)
    
    def analyze_all_events(self):
        
        e = self.e
        self.q_induced = []
        self.q_direct = []
        start_time = time.time()
        bar_length = 50.
        n = len(self.split_data)
        for i,event in enumerate(self.split_data):
            e.set_track(event[:,[2,3,4]], event[:,5] * 1e6)
            xyz, polarity = e.get_charge(fQ=self.fQ)
            qi = e.get_induced_charge_on_qubits(xyz, polarity)
            # qi = e.get_induced_charge_on_qubits(np.full((1,1),(0,0,0)), [1])
            qd = e.get_charge_on_qubits(xyz, polarity)
            self.q_induced.append(qi)
            self.q_direct.append(qd)
            
            # Update progress bar
            progress = (i+1.)/n
            finish_time = time.localtime(
                start_time + n * (time.time() - start_time) / (i+1) )
            bartxt = '[' + '#'*int(np.floor(bar_length*progress)) \
                        + ' '*int(np.ceil(bar_length*(1.-progress))) \
                        + '] {}% '.format(int(100.*progress)) \
                        + time.strftime("%H:%M:%S %m/%d/%Y\r", finish_time)
            sys.stdout.flush()
            sys.stdout.write(bartxt)
            sys.stdout.flush()
            
        self.q_induced = np.array(self.q_induced).reshape(-1,4)
        self.q_direct = np.array(self.q_direct).reshape(-1,4,2)
    
    def plot_qq(self, q1, q2, ax=None, q_induced=None):
        """ Shows 2d plot showing the aliased charge induced from all events, with
        one axis as q1 and the other axis q2. """
    
        if q_induced is None:
            q_induced = self.q_induced
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.plot( noiselib.alias(self.q_induced[:,q1-1]),
                 noiselib.alias(self.q_induced[:,q2-1]), '.' )
        ax.set_xlabel('Q{} aliased charge [e]'.format(q1))
        ax.set_ylabel('Q{} aliased charge [e]'.format(q2))
        ax.set_title('Q{} - Q{}'.format(q1,q2))
        ax.set_aspect(1)
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,0.5)
        plt.draw()
        plt.pause(0.05)
    
    def plot_charge_hist(self, q):
        """ Shows histogram of the induced charge on qubit q for ALL events. """
        
        fig, ax = plt.subplots(1,1)
        ax.hist( noiselib.alias(self.q_induced[:,q-1]), bins=50, range=(-0.5,0.5) )
        ax.set_xlabel('Q{} aliased charge [e]'.format(q))
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of induced charge from all events on Q{}'.format(q))
        plt.draw()
        plt.pause(0.05)
    
    def plot_charge_hist_from_impact(self, q):
        """ Show histogram of charge induced on qubit q from CURRENT event. """
    
        event = self.events[ self.i % len(self.split_data) ]
        charge_xyz, charge_polarity = event.get_charge(fQ=self.fQ)
        qs = event.get_induced_charge_on_qubits(charge_xyz, charge_polarity, sum=False)
        qs = qs[q-1,:]
        fig, ax = plt.subplots(1,1)
        ax.hist([qs[charge_polarity<0], qs[charge_polarity>0]], bins=50)
        ax.set_xlabel('Q{} unaliased charge [e]'.format(q))
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of induced charge from event {} on Q{}'.format( 
                            self.i % len(self.split_data), q ))
        ax.legend(['electrons','holes'])
        plt.draw()
        plt.pause(0.05)
    
    def count(self, q1, q2, quadrant, thresh=0.1, q_induced=None):
    
        if q_induced is None:
            q_induced = self.q_induced
        
        e1 = noiselib.alias(q_induced[:,q1-1])
        e2 = noiselib.alias(q_induced[:,q2-1])
        
        if quadrant == 1:
            return np.sum( (e1 > thresh) & (e2 > thresh) )
        elif quadrant == 2:
            return np.sum( (e1 < -thresh) & (e2 > thresh) )
        elif quadrant == 3:
            return np.sum( (e1 < -thresh) & (e2 < -thresh) )
        elif quadrant == 4:
            return np.sum( (e1 > thresh) & (e2 < -thresh) )
    
    def get_correlation(self, q1, q2, thresh=0.1, q_induced=None):
    
        if q_induced is None:
            q_induced = self.q_induced
        
        e1 = noiselib.alias(q_induced[:,q1-1])
        e2 = noiselib.alias(q_induced[:,q2-1])
        return 1.*np.sum( (np.abs(e1) > thresh) & (np.abs(e2) > thresh) ) / \
                  np.sum( (np.abs(e1) > thresh) | (np.abs(e2) > thresh) )
        

if __name__ == '__main__':
    app = Controller(sys.argv)
    if hasattr(app, 'view'):
        sys.exit(app.exec_())