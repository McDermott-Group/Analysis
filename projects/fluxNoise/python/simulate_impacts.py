import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
import matplotlib.pyplot as plt
import noiselib
from numba import jit


class ImpactEvent(object):
    
    def __init__(self, pdfs, path, energy, qubit_xys, charge_map):
        self.PDFs = pdfs
        self.path = path
        self.energy = energy
        self.qubit_xys = qubit_xys
        self.charge_map = charge_map
    
    def set_path(self, path, energy):
        self.path = path
        self.energy = energy
        if hasattr(self, 'charge'):
            del self.charge
    
    def get_path(self):
        c = self.energy
        colorFloats = np.interp(c, (c.min(), c.max()), (0,1))
        color = cm.autumn(colorFloats)
        return self.path, color
    
    def get_charge(self, fQ=1., cache=False):
        if hasattr(self, 'charge'):
            return self.charge
        xyz = []
        polarity = []
        for i in range(len(self.path)):
            x,y,z = self.path[i]
            res = self.diffuseCharge(x,y,z,self.energy[i],fQ=fQ,epseh=3.6,verbose=False)
            polarity += [np.full(res['electrons'].shape[1], -1.), 
                         np.full(res['holes'].shape[1], 1.) ]
            xyz += [res['electrons'].T, res['holes'].T]
        xyz = np.concatenate(xyz)
        polarity = np.concatenate(polarity)
        if cache:
            self.charge = xyz, polarity
            print xyz.nbytes/1e6
        return xyz, polarity
    
    def get_charge_on_qubits(self, charge_xyz, charge_polarity):
        q = []
        for i,xy_q in enumerate(self.qubit_xys):
            dx,dy,dz = (charge_xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
            r = np.hypot(dx,dy)
            qs = (r < 0.07) & (dz == 0)
            q.append( np.sum(qs*charge_polarity) )
        return q
    
    def get_induced_charge_on_qubits(self, charge_xyz, charge_polarity):
        q = []
        for i,xy_q in enumerate(self.qubit_xys):
            dx,dy,dz = (charge_xyz - (xy_q[0], xy_q[1], 0.375/2.)).T
            r = np.hypot(dx,dy)
            qs = self.charge_map(r,-dz)
            # if i == 0:
                # # plt.hist([qs[charge_polarity<0], qs[charge_polarity>0]],bins=100)
                # # plt.draw()
                # # plt.pause(0.05)
                # plt.figure()
                # plt.hist2d(charge_xyz[:,0], charge_xyz[:,1], bins=100)
                # plt.pause(0.05)
                # plt.figure()
            q.append( np.sum(qs*charge_polarity) )
        return q
    
    def potential(self, r, zr, p, zp, z_shift):
        """Calculate the potential at a point (p,0,zp) due to a 
        ring of charge Q=1 centered at (0,0,zr) with radius r."""
        z = zr - zp + z_shift
        rd = r - p
        m = -4.*r*p/(z**2+rd**2)
        ellk = 1./np.sqrt(1.-m)*ellipk(-m/(1.-m)) # Imaginary-modulus transformation
        v = 1./(2*np.pi) * 4./np.sqrt(z**2.+rd**2.)*ellk
        return v
    
    def generate_induced_charge_map(self):
        r_inner, r_outer = 70, 90.5
        r_grid, z_grid = np.arange(1,550,80), np.arange(0,375,100)
        r_rings = np.concatenate([np.arange(1,r_inner), 
                                  np.arange(r_outer,500,1), [0]])
                                  
        z_rings = np.zeros(r_rings.size)
        v_rings = np.concatenate([ np.zeros(r_rings.size-1), [1] ])
        q_rings = np.zeros(r_rings.size)
        C_inv = np.zeros( (r_rings.size,r_rings.size) )
        for i in range(r_rings.size-1):
            C_inv[:,i] = self.potential(r_rings[i], z_rings[i], r_rings, z_rings, 0.1)
        induced_charge = np.zeros( (r_grid.size, z_grid.size) )
        for i,z in enumerate(z_grid):
            for j,r in enumerate(r_grid):
                print i,j
                r_rings[-1] = r
                z_rings[-1] = z
                C_inv[-1,:] = C_inv[:,-1] = self.potential(r_rings[-1], z_rings[-1], r_rings, z_rings, 0.1)
                # q_rings,_,_,_ = np.linalg.lstsq(C_inv, v_rings)
                q_rings = np.dot( np.linalg.inv(C_inv), v_rings )
                induced_charge[j,i] = ( np.sum(q_rings[:np.arange(1,r_inner).size]) /
                                        np.sum(q_rings) )
        return induced_charge
    
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
            zf = z+sim_distz
            
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
    
    def __init__(self):
        super(MainWindow, self).__init__()
        # self.setGeometry(300, 300, 250, 150)
        self.show()
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
        
        self.qubit_xys = [( 1.564,  0.570),
                          ( 1.564, -0.070),
                          (-1.564, -0.080),
                          (-1.564, -0.420)]
        for x,y in self.qubit_xys:
            self.draw_qubit(x,y)
        
        PDFs = np.load('sim_data/ChargePDFs.npy',allow_pickle=True)
        PDFs = PDFs.tolist()
        
        q_map = noiselib.loadmat('sim_data/charge_map.mat')
        q = q_map['charge_mat550_5000']
        q[~np.isfinite(q)] = -1. # right below qubit
        # charge_map = interp2d( q_map['r']/1000., q_map['z']/1000., q, copy=False)
        # charge_map = np.vectorize(charge_map) # much slower
        def map_charge(r, z, q, rp, zp): 
            # my own bilinear interpolation, much faster
            # https://math.stackexchange.com/questions/3230376/interpolate-between-4-points-on-a-2d-plane
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
        charge_map = lambda rp,zp: map_charge( q_map['r']/1000., q_map['z']/1000., q, rp, zp)
        
        data = np.loadtxt('sim_data/Muons.txt', skiprows=1)
        # data = np.loadtxt('sim_data/Gamma.txt', skiprows=1)
        self.split_data = np.split(data, np.where(np.diff(data[:,0]))[0]+1)
        self.events = []
        for event in self.split_data:
            e = ImpactEvent(PDFs, event[:,[2,3,4]], event[:,5] * 1e6,
                            self.qubit_xys, charge_map)
            self.events.append(e)
        
        self.i = 0
        self.plot_trace(0)
        
        # e = ImpactEvent(PDFs, 1, 1, self.qubit_xys, charge_map)
        # ps = []
        # for event in self.split_data:
            # e.set_path(event[:,[2,3,4]], event[:,5] * 1e6)
            # xyz, polarity = e.get_charge()
            # p = e.get_induced_charge_on_qubits(xyz, polarity)
            # print p
            # ps.append(p)
        # ps = np.array(ps).reshape(-1,4)
        # plt.hist(ps.T)
        # plt.draw()
        # plt.pause(0.05)
        
        # threads = []
        # for e in self.events:
            # thread = threading.Thread(target=e.get_charge)
            # threads.append(thread)
            # thread.start()
        # for thread in threads:
            # thread.join()
    
    def chip(self):
        h = 0.375/2.
        L = 6.25/2.
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
    
    def draw_chip_lines(self):
        h = 0.375/2.
        L = 6.25/2.
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
    
    def draw_qubit(self, x, y):
        r_outer = 0.0905
        r_inner = 0.07
        h = 0.375/2.
        
        outer_circle = self._circle((x,y), r_outer)
        outer_circle = np.concatenate([outer_circle, [outer_circle[0]]])
        self.addItem( gl.GLLinePlotItem(pos=outer_circle) )
        
        inner_circle = self._circle((x,y), r_inner)
        vertexes = np.concatenate([[[x,y,h]], inner_circle])
        pts = range(1,len(vertexes))
        faces = np.array([ [0,pts[i],pts[i+1]] for i in range(-1,len(pts)-1) ])
        self.addItem( gl.GLMeshItem(vertexes=vertexes, faces=faces) )
    
    def _circle(self, xy, r):
        """Returns points for a circle."""
        n = 100;
        h = 0.375/2.
        x = xy[0] + r*np.cos(np.linspace(0,2*np.pi,n))
        y = xy[1] + r*np.sin(np.linspace(0,2*np.pi,n))
        h = np.full(n, h)
        return np.transpose([x,y,h])
        
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Left:
            self.i -= 1
        elif event.key() == QtCore.Qt.Key_Right:
            self.i += 1
        else:
            return
        self.plot_trace(self.i%len(self.split_data))
    
    def plot_trace(self, i):
        event = self.events[i]
        pts, color = event.get_path()
        self.track.setData(pos=pts,
                           color=color)
        pts, polarity = event.get_charge()
        self.charge.setData(pos=pts,
                            color=cm.PiYG(polarity))
        print event.get_induced_charge_on_qubits(pts, polarity)
        print event.get_charge_on_qubits(pts, polarity)


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    w = MainWindow()
    sys.exit(app.exec_())


# PDFs = np.load('sim_data/ChargePDFs.npy',allow_pickle=True).tolist()
# PDFs = PDFs.tolist()
# q_map = noiselib.loadmat('sim_data/charge_map.mat')
# r, z = q_map['r']/1000., q_map['z']/1000.
# q = q_map['charge_mat550_5000']
# charge_map = interp2d( r, z, q, copy=False)
# charge_map = np.vectorize(charge_map)
# data = np.loadtxt('sim_data/Muons.txt', skiprows=1)
# # data = np.loadtxt('sim_data/Gamma.txt', skiprows=1)
# split_data = np.split(data, np.where(np.diff(data[:,0]))[0]+1)

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import simulate_impacts
# reload(simulate_impacts)
# from simulate_impacts import ImpactEvent
# PDFs = np.load('sim_data/ChargePDFs.npy',allow_pickle=True).tolist()
# e = ImpactEvent(PDFs, np.array((0,0,0)).reshape(1,3), np.array((.1e6,)), 0, 1)
# xyz, pol = e.get_charge()
# x,y,z = xyz.T
# fig, (ax1,ax2) = plt.subplots(1,2)
# h1 = ax1.hist2d(x[pol<0], y[pol<0], bins=100, norm=mpl.colors.LogNorm());
# h2 = ax2.hist2d(x[pol>0], y[pol>0], bins=100, norm=mpl.colors.LogNorm());
# ax1.set_aspect(1); ax2.set_aspect(1)
# ax2.get_yaxis().set_visible(False)
# ax1.set_title('electrons'); ax2.set_title('holes');
# fig.colorbar(h1[3], ax=(ax1,ax2)); #plt.colorbar();
# fig.suptitle('XY'); plt.draw(); plt.pause(0.05)
