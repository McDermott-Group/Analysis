try:
    import OpenGL as ogl
    try:
        import OpenGL.GL   # this fails in <=2020 versions of Python on OS X 11.x
    except ImportError:
        print('Drat, patching for Big Sur')
        from ctypes import util
        orig_util_find_library = util.find_library
        def new_util_find_library( name ):
            res = orig_util_find_library( name )
            if res: return res
            return '/System/Library/Frameworks/'+name+'.framework/'+name
        util.find_library = new_util_find_library
except ImportError:
    pass

from importlib import reload
import numpy as np
import matplotlib.cm as cm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
from scipy import constants
import matplotlib.pyplot as plt
import noiselib
# from numba import jit
import time
# import multiprocessing as mp
# import pathos.pools as pp

# base_path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
base_path = '/Volumes/smb/mcdermott-group/data/fluxNoise2/sim_data'

class ImpactEvent(object):

    def __init__(self, pdfs, charge_map, qubit_xys, track=None, energy=None):

        if len(pdfs) > 1:
            self.PDFs = {'electrons': pdfs[0],
                         'holes': pdfs[1] }
        else:
            self.PDFs = {'electrons': pdfs[0],
                         'holes': pdfs[0] }
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

    def get_charge(self, fQ=1., cache=False, concatenate=True):
        """ Generates and returns positions and polarities of charges from track
        and energy. """

        if hasattr(self, 'charge'):
            return self.charge
        xyz = []
        polarity = []
        for i in range(len(self.track)):
            x,y,z = self.track[i]
            res = self.diffuseCharge(x,y,z,self.energy[i],fQ=fQ,epseh=3.6,verbose=False)
            polarity += [np.concatenate([np.full(res['electrons'].shape[1], -1.),
                                        np.full(res['holes'].shape[1], 1.) ])]
            xyz += [np.concatenate([res['electrons'].T, res['holes'].T])]
        if concatenate:
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

    def get_induced_qubit_rotation(self, charge_xyz, sum=True):
        """ Calculate the induced rotation on each qubit from charges in the method
        get_induced_charge_on_qubits() above.  If sum==True,
        the induced rotations will be summed before being returned, giving the total
        rotation.  If False, an array of induced rotations will be returned
        for each qubit. """

        qs = self.get_induced_charge_on_qubits(charge_xyz, 1, sum=False)
        qs = noiselib.alias(qs, 1.)
        Ec = 2 * np.pi * 235e6
        w_01 = 2 * np.pi * 4.1e9
        if sum:
            qs2 = np.sum(qs**2, axis=1)
        else:
            qs2 = qs**2
        return 2*np.sqrt(Ec/w_01) * np.sqrt(qs2)

    def getPDF(self, z, charge_type='electrons'):

        zcuts = self.PDFs[charge_type]['zarray']
        zind = np.argmin(np.abs(z-zcuts))
        label = (charge_type,zcuts[zind],charge_type)
        if label not in self.PDFs:
            Hn = self.PDFs[charge_type][zcuts[zind]][charge_type]
            Hr = Hn.ravel()
            self.PDFs[label] = np.cumsum(Hr)/np.sum(Hr), Hn.shape
        return self.PDFs[label]

    def diffuseCharge(self, x, y, z, energy, fQ = 1.0, epseh=3.6, verbose=True):

        PDFs = self.PDFs['electrons']

        neh = int(np.floor(energy/epseh)*fQ)
        if(verbose):
            print('Event Energy: '+str(energy)+' eV')
            print('Electron-Hole Pairs: '+str(neh))

        results=dict()
        for charge_type in ['electrons','holes']:
            # Hn=self.getPDF(z,charge_type=charge_type)
            # Hr = Hn.ravel()
            # normcumsum = np.cumsum(Hr)/np.sum(Hr)
            # vals = np.random.choice(Hr.size, p=Hr.astype('float64'), size=neh)
            normcumsum, shape = self.getPDF(z,charge_type=charge_type)
            r = np.random.random(neh)
            vals = np.searchsorted(normcumsum, r)
            inds = np.unravel_index(vals, shape)

            sim_distx = PDFs['x'][inds[0]]
            sim_disty = PDFs['y'][inds[1]]
            sim_distz = PDFs['z'][inds[2]]

            xf = x+sim_distx
            yf = y+sim_disty
            zf = sim_distz

            xf += +np.random.rand(len(xf))*PDFs['dx'] - PDFs['dx']/2.0
            yf += +np.random.rand(len(yf))*PDFs['dy'] - PDFs['dy']/2.0
            zf += +np.random.rand(len(zf))*PDFs['dz'] - PDFs['dz']/2.0

            zmax = PDFs['thickness']/2.
            zmin = -PDFs['thickness']/2.
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
        self.setGeometry(300, 300, 1600*3/4, 900*3/4)
        self.show()
        self.setBackgroundColor('w')
        self.model = model

        # g = gl.GLGridItem(color=pg.glColor('w'))
        # g.setSize(6.25,6.25,0.375)
        # self.addItem(g)

        self.track = gl.GLScatterPlotItem( size=7 )
        self.track.setGLOptions('translucent')
        self.addItem(self.track)

        self.charge = gl.GLScatterPlotItem( size=3 )
        self.charge.setGLOptions('translucent')
        self.addItem(self.charge)

        # ax = gl.GLAxisItem()
        # ax.setSize(10,10,10)
        # self.addItem(ax)

        # self.setBackgroundColor('w')
        self.draw_chip_lines()
        # self.draw_chip_faces()
        self.draw_chip()

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
        self.addItem( gl.GLMeshItem(vertexes=vertexes, faces=faces) )

    def draw_chip(self):
        img_dat = plt.imread('die_stitched.png')
        rbga = np.dstack([img_dat, np.full((1089,1089),0.7)])
        img = gl.GLImageItem(rbga)
        img.translate(-6.25/2, -6.25/2, 0)
        # img.rotate(90, 0,0,1)
        # img.rotate(-90, 0,1,0)
        self.addItem(img)
        print('done')

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
        self.addItem( gl.GLLinePlotItem(pos=vertexes, color=pg.glColor('k'), width=2, antialias=True) )

    def draw_qubit(self, x, y, r_inner=0.07, r_outer=0.0905, h=0.375/2.):

        outer_circle = self._circle((x,y), r_outer)
        outer_circle = np.concatenate([outer_circle, [outer_circle[0]]])
        self.addItem( gl.GLLinePlotItem(pos=outer_circle, color=pg.glColor('k'), width=2, antialias=True))

        inner_circle = self._circle((x,y), r_inner)
        vertexes = np.concatenate([[[x,y,h]], inner_circle])
        pts = range(1,len(vertexes))
        faces = np.array([ [0,pts[i],pts[i+1]] for i in range(-1,len(pts)-1) ])
        self.addItem( gl.GLMeshItem(vertexes=vertexes, faces=faces, color=pg.glColor('k')) )

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
        elif event.key() == 83:
            if hasattr(self.model, '_save_frame'):
                self.model._save_frame()
        else:
            return

    def plot_track(self, track_xyz, track_color, charge_xyz, charge_polarity):

        self.track.setData(pos=track_xyz,
                           color=track_color)
        self.charge.setData(pos=charge_xyz,
                            color=cm.bwr(charge_polarity)) # PiYG


class Controller(QtGui.QApplication):

    def __init__(self, sys_argv, plot=True, calc_all=False,
                 pdfs_file='{}/ChargePDFs_{}.npy'.format(base_path,100),
                 event_files=[base_path+'/Gamma.txt',
                              base_path+'/Gamma_10deg.txt',
                              base_path+'/Gamma_10deg_pt2.txt'],
                 fQ=1.):

        super(Controller, self).__init__(sys_argv)

        self.qubit_xys = [( 1.564,  0.570),
                          ( 1.564, -0.070),
                          (-1.564, -0.080),
                          (-1.564, -0.420)]

        self.fQ = fQ
        self.event_files = event_files
        self.txt_data = []

        if type(pdfs_file) != list:
            pdfs_file = [pdfs_file]
        PDFs = [np.load(p_file, allow_pickle=True).tolist() for p_file in pdfs_file]

        charge_map = self.load_charge_map(base_path+'/charge_map2.mat')

        self.event = ImpactEvent(PDFs, charge_map, self.qubit_xys)

        self.i = -1
        self.n_events = self.get_num_events()

        if plot:
            self.view = MainWindow(self, self.qubit_xys)
            self.show_next()

        if calc_all:
            self.e = ImpactEvent(PDFs, 1, 1, self.qubit_xys, charge_map)
            self.thread = threading.Thread(target=self.analyze_all_events)
            self.thread.start()
            # self.analyze_all_events()

    def load_charge_map(self, path):

        q_map = noiselib.loadmat(path)
        q = q_map['charge_mat']
        q[~np.isfinite(q)] = -1. # right below center of qubit
        charge_map = lambda rp,zp: noiselib.interp2d( q_map['r']/1000.,
                                                      q_map['z']/1000.,
                                                      q, rp, zp )
        return charge_map

    def load_track_data(self, i):

        txt_data = []
        last_index = 0
        index = 0
        for path in self.event_files:
            with open(path, 'r') as f:
                for j,line in enumerate(f):
                    if j==0:
                        continue
                    new_index = int(line.split()[0])
                    if new_index > last_index:
                        index += 1
                    if index == i:
                        txt_data.append(line)
                    elif index > i:
                        data = np.loadtxt(txt_data)
                        if len(data.shape) == 1:
                            data = data.reshape(1,data.size)
                        return data
                    last_index = new_index

    # def load_track_data(self, i):

        # txt_data = self.txt_data
        # data = []
        # for path in self.event_files:
            # data.append( np.loadtxt(path, skiprows=1) )
        # data = np.concatenate(data)
        # split_data = np.split(data, np.where(np.diff(data[:,0]))[0]+1)
        # return split_data[i]

    def get_next_track(self):

        txt_data = []
        last_index = 0
        for path in self.event_files:
            with open(path, 'r') as f:
                for j,line in enumerate(f):
                    if j==0:
                        continue
                    new_index = int(line.split()[0])
                    if new_index > last_index:
                        data = np.loadtxt(txt_data)
                        if len(data.shape) == 1:
                            data = data.reshape(1,data.size)
                        txt_data = [line]
                        yield data
                    else:
                        txt_data.append(line)
                    last_index = new_index

    def get_num_events(self):

        index = 0
        for path in self.event_files:
            f_index = 0
            with open(path, 'r') as f:
                for j,line in enumerate(f):
                    if j==0:
                        continue
                    f_index = int(line.split()[0])
            index += f_index
        return index + 1

    def show_previous(self):

        self.i -= 1
        self.show_event()

    def show_next(self):

        self.i += 1
        self.show_event()

    def sim_event(self):

        e = self.load_track_data( self.i % self.n_events )
        self.event.set_track(e[:,[2,3,4]], e[:,5] * 1e6)
        track_xyz, track_color = self.event.get_track()
        charge_xyz, charge_polarity = self.event.get_charge(fQ=self.fQ)
        return track_xyz, track_color, charge_xyz, charge_polarity

    def show_event(self):
        track_xyz, track_color, charge_xyz, charge_polarity = self.sim_event()
        self.view.plot_track(track_xyz, track_color, charge_xyz, charge_polarity)

    def analyze_all_events(self):

        self.q_induced = []
        self.q_direct = []
        self.rot_induced = []
        start_time = time.time()
        bar_length = 50.
        track_gen = self.get_next_track()
        n = self.n_events
        for i in range(n):
            try:
                e = track_gen.next()
            except StopIteration:
                print('extra event attempt')
                continue
            self.event.set_track(e[:,[2,3,4]], e[:,5] * 1e6)
            xyz, polarity = self.event.get_charge(fQ=self.fQ)
            qi = self.event.get_induced_charge_on_qubits(xyz, polarity)
            roti = self.event.get_induced_qubit_rotation(xyz)
            # qi = self.event.get_induced_charge_on_qubits(np.full((1,1),(0,0,0)), [1])
            qd = self.event.get_charge_on_qubits(xyz, polarity)
            self.q_induced.append(qi)
            self.rot_induced.append(roti)
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
        self.rot_induced = np.array(self.rot_induced).reshape(-1,4)
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

        e = self.load_track_data( self.i % self.n_events )
        self.event.set_track(e[:,[2,3,4]], e[:,5] * 1e6)
        charge_xyz, charge_polarity = self.event.get_charge(fQ=self.fQ)
        qs = self.event.get_induced_charge_on_qubits(charge_xyz, charge_polarity,
                                                     sum=False)
        qs = qs[q-1,:]
        fig, ax = plt.subplots(1,1)
        ax.hist([qs[charge_polarity<0], qs[charge_polarity>0]], bins=50)
        ax.set_xlabel('Q{} unaliased charge [e]'.format(q))
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of induced charge from event {} on Q{}'.format(
                            self.i % self.n_events, q ))
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
        both_jump = np.sum( (np.abs(e1) > thresh) & (np.abs(e2) > thresh) )
        either_jump = np.sum( (np.abs(e1) > thresh) | (np.abs(e2) > thresh) )
        try:
            corr = 1.*both_jump / either_jump
        except ZeroDivisionError:
            corr = np.nan
        try:
            dcorr = 1.*both_jump / either_jump * np.sqrt(1./both_jump + 1./either_jump)
        except ZeroDivisionError:
            dcorr = np.nan
        return corr, dcorr



if __name__ == '__main__':
    app = Controller(sys.argv)
    if hasattr(app, 'view'):
        sys.exit(app.exec_())
