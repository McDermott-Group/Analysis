import numpy as np
import matplotlib.cm as cm
#from pyqtgraph.Qt import QtCore, QtGui
#import pyqtgraph as pg
#import pyqtgraph.opengl as gl
import sys, os
import threading
from scipy.special import ellipk
from scipy.interpolate import interp2d, griddata
from scipy import constants
import matplotlib.pyplot as plt
import noiselib
from importlib import reload
# from numba import jit
import time
import impact_lib
reload(impact_lib)
from impact_lib import *
import pickle
import gc
import itertools

#path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data'
path = '/Volumes/smb/mcdermott-group/data/fluxNoise2/sim_data'
thresh = 0.1
hit_type = 'muons'
L_list = [100,200,300,400,500,600,700,800,900,1000]
fq_list = [ 1.,0.5,0.2,0.1 ]
L, fq = 300, 0.2


def load_pdfs(L):
    # load charge pdf
    if type(L) in (list,tuple):
        pdfs_file = ['{}/ChargePDFs_{}r_vhs_fine.npy'.format(path,l) for l in L]
    else:
        pdfs_file = '{}/ChargePDFs_{}r_vhs_fine.npy'.format(path,L)
    return pdfs_file

def load_hits(hit_type):
    # load generated hits
    if hit_type == 'gammas':
        efiles = [path+'/Gamma.txt',
                  path+'/Gamma_10deg.txt',
                  path+'/Gamma_10deg_pt2.txt']
    elif hit_type == 'muons':
        efiles = [path+'/Muons.txt']#, path+'/Muons2.txt']
    return efiles


pdfs_file = load_pdfs(L)
efiles = load_hits(hit_type)


class AnimationController(Controller):

    def __init__(self, sys_argv, plot=False, calc_all=False,
                 pdfs_file='{}/ChargePDFs_{}.npy'.format(base_path,100),
                 event_files=[base_path+'/Gamma.txt',
                              base_path+'/Gamma_10deg.txt',
                              base_path+'/Gamma_10deg_pt2.txt'],
                 fQ=1.):

        super(AnimationController, self).__init__(sys_argv, fQ=fQ, calc_all=False,
                          plot=False, pdfs_file=pdfs_file, event_files=efiles)

        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        self.view = MainWindow(self, self.qubit_xys)

        self.v_track = 0.2
        self.v_charge = 0.05

        self.i = 78 # 8, 19, 78
        self.animate_event()

    def show_event(self):
        if hasattr(self,'tt'):
            self.tt.stop()
        self.animate_event()

    def sim_event(self):

        print(self.i)
        e = self.load_track_data( self.i % self.n_events )
        self.event.set_track(e[:,[2,3,4]], e[:,5] * 1e6)
        track_xyz, track_color = self.event.get_track()
        charge_xyz, charge_polarity = self.event.get_charge(fQ=self.fQ, concatenate=False)
        return track_xyz, track_color, charge_xyz, charge_polarity

    def animate_event(self):
        track_xyz, track_color, charge_xyz, charge_polarity = self.sim_event()

        ncharge = [charge.shape[0] for charge in charge_xyz]
        indcs = np.cumsum(ncharge)
        indcs_ext = np.concatenate([np.full((n,3),i) for i,n in enumerate(ncharge)])
        track_xyz_ext = np.concatenate([np.full((n,3),track_xyz[i]) for i,n in enumerate(ncharge)])
        charge_xyz = np.concatenate(charge_xyz)
        charge_polarity = np.concatenate(charge_polarity)
        charge_xyz_dif = charge_xyz - track_xyz_ext
        d_track = np.concatenate([(0,0,0), np.cumsum(np.linalg.norm(np.diff(track_xyz,axis=0),axis=1)) ])
        d_track_ext = np.concatenate([np.full((n,3),d_track[i]) for i,n in enumerate(ncharge)])
        self.track_charge_data = (track_xyz, track_color,
                                  charge_xyz, charge_polarity,
                                  indcs, indcs_ext,
                                  track_xyz_ext, charge_xyz_dif,
                                  d_track, d_track_ext)
        self.view.orbit(0,-0)

        self.t = 20. # 1.
        self.img_i = 0
        self.update_animation()
        if False: # animate
            self.tt = QtCore.QTimer()
            self.tt.timeout.connect(self.update_animation)
            self.tt.start(50)

        # self.view.plot_track(track_xyz, track_color,
        #                      charge_xyz, charge_polarity)

    def update_animation(self):
        self.t += 0.05
        (track_xyz, track_color, charge_xyz, charge_polarity,
        indcs, indcs_ext, track_xyz_ext, charge_xyz_dif, d_track, d_track_ext) = self.track_charge_data
        # i = min( int(np.floor(self.v_track*self.t)), track_xyz.shape[0]-1 )
        j = min( np.searchsorted(d_track,self.v_track*self.t), track_xyz.shape[0]-1 )
        i = j
        # mag = np.minimum( np.abs(i-indcs_ext)/self.v_track * self.v_charge * np.abs(charge_xyz_dif), np.abs(charge_xyz_dif) )
        mag = np.clip( (self.t - d_track_ext/self.v_track-10) * self.v_charge * np.abs(charge_xyz_dif), 0., np.abs(charge_xyz_dif) )
        charge_xyz_t = track_xyz_ext + mag * np.sign(charge_xyz_dif)
        self.view.plot_track(track_xyz[0:i], track_color[0:i],
                             charge_xyz_t[0:indcs[i]], charge_polarity[0:indcs[i]])
        self.view.orbit(-0.1,-0.025)
        d = self.view.cameraPosition().length()
        self.view.setCameraPosition(distance=d-0.005)
        # self._save_frame()
        if np.all( (self.t - d_track_ext/self.v_track) * self.v_charge * np.abs(charge_xyz_dif) >= np.abs(charge_xyz_dif) ):
            self.tt.stop()
        # print(self.view.opts['elevation'])
        # print(self.view.opts['azimuth'])

    def _save_frame(self):
        # self.view.grabFrameBuffer().save('animation/frame_{:03d}.png'.format(self.img_i))
        frame = self.view.grabFrameBuffer()
        print('saving frame')
        print(frame)
        frame.save('stills/frame_{:03d}_{:03d}.png'.format(self.img_i, self.i))
        self.img_i += 1



if __name__ == '__main__':
    app = AnimationController(sys.argv, fQ=fq, calc_all=False, plot=True,
                     pdfs_file=pdfs_file, event_files=efiles)
    sys.exit(app.exec_())
