import numpy as np

class QP0D(object):
    """
    This is to analyze the 0D QP Dynamics
    No spatial
    Just Energy
    Follow Martinis/Leanader/McDermott
    """

    def __init__(self):
        ### some constants
        self.Delta = 1  # Energy gap
        self.tau0 = 438  # e-p coupling time, units of ns, use Kaplan's value
        self.Dn = 6e-2  # units of um^2/ns -- use 6 here for aluminum


        ###
        self.n = None    #
        self.e_grid = None   # from

    def import_data(self):
        self._simulationSetup()
        self._generateEneryGrid()
        self._qpOccupationProbability()

    def _simulationSetup(self):
        self.nx = 1     # points in space grid
        self.ne = 100   # points in energy grid

        self.time_tot = 0.001 * 1e-3    # total time for simulation
        self.step = 1e-4  # sets time step precision
        self.dt = self.tau0 * self.step # each step's time
        self.nt = (self.time_tot*1e9/self.dt)   # points in time grid

    def _generateEneryGrid(self):
        e_max = 10.0 # units of Delta, upper limit of energy space
        e_grid = np.logspace(0, 4, num=self.ne, base=2)
        e_inj = 5   # injection energy at 5 Delta
        idx_e_upper = (np.abs(e_grid - e_inj)).argmin()
        e_grid = e_grid[:idx_e_upper + 5]  # useful energy grid
        self.e_grid = e_grid
        self.n = np.zeros(len(e_grid))

    def _qpOccupationProbability(self):
        e_grid = self.e_grid

        Gamma = 1e-3  # Dynes broadening
        broadened = complex(0, -Gamma)
        Dynes = e_grid + broadened  # to avoid divergence at the gap
        rho = Dynes / np.sqrt(Dynes ** 2 - Delta ** 2)
        rho = rho.real
        self.rho = rho

    def _generateScatteringMatrix(self):
        ne = self.ne
        rho = self.rho
        G_s = np.zeros((ne, ne))

        ### kTc = \Delta/1.76
        for i in range(ne):
            for j in range(ne):
                # this is scattering from j to i
                G_s[i, j] = s * (1.76 ** 3.0) * ((e[i] - e[j]) ** 2.0) * (1.0 - Delta ** 2.0 / (e[i] * e[j])) \
                          * np.heaviside(e[j] - e[i], 0.5) * rho[j]
                # np.heaviside(e[j] - e[i], 0.5) the phonon occupation factor = 0, 0.5, 1.0 in this case

                # print('PhononTemp=', np.heaviside(e[j] - e[i], 0.5))

