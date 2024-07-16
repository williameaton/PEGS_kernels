from abc import ABC, abstractmethod
from wetools import obspy_gen_mpl
import numpy as np
import obspy

class PegsBase(ABC):

    @abstractmethod
    def __init__(self, filepath, **kwargs):

        # Default kwargs:
        defaultKwargs = {'stations': ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO'],
                         'cutoff'  : {'MAJO': 57.95,  'INU':  77.77, 'MDJ':  163.81,  'INCN': 177.7 ,  'NE93': 237.43, 'BJT':  276.4 , 'MA2':  292.3 ,'XAN':  343.65,  'ULN': 344. },
                         'scale'   : 1e9}
        kwargs = defaultKwargs | kwargs

        # Store inputs
        self.fpath = filepath

        # Possible kwargs:
        self.stns  = kwargs['stations']
        self.scale = kwargs['scale']
        self.cutoff = kwargs['cutoff']

        # Setup
        self.pegs = {}
        self.grav = {}
        self.zacc = {}
        self.time = {}
        self.real = {'time' : {},
                     'pegs' : {}}

        self.traces = {'time': self.time,
                       'pegs': self.pegs,
                       'grav': self.grav,
                       'zacc': self.zacc}

        self.z_dir = {'grav': 1,
                      'zacc': 1}

        self.nstns = len(self.stns)




    def load(self):
        self._load_data('grav')
        self._load_data('zacc')
        self._compute_pegs()


    def _load_data(self, dtype):

        for stn in self.stns:
            # Load to numpy arrays
            t, y = obspy_gen_mpl(obspy.read(f'{self.fpath}/{dtype}/{stn}_proc.sac')[0])

            # Add to relevant array:
            self.traces['time'][stn] = np.array(t)
            self.traces[dtype][stn] = np.array(y*self.z_dir[dtype])* self.scale


    def _compute_pegs(self):
        for stn in self.stns:
            self.pegs[stn] = self.grav[stn] - self.zacc[stn]


    def set_z_dir(self, dtype, value):
        if np.abs(value)!=1:
            raise ValueError('value must be +-1')
        if dtype!="grav" and dtype!="zacc":
            raise ValueError(f"dtype must be 'grav' or 'zacc'")

        self.z_dir[dtype] = value

    def cut_to_plot(self, plusval=35):
        for stn in self.stns:
            t = self.time[stn]
            m = t < self.cutoff[stn] + plusval

            self.time[stn] = t[m]

            self.grav[stn] = self.grav[stn][m]
            self.zacc[stn] = self.zacc[stn][m]


    def load_real_data(self, fpath):

        for stn in self.stns:
            t, y = obspy_gen_mpl(obspy.read(f'{fpath}/{stn}_acc_Z')[0])

            # Add to relevant array:
            self.real['time'][stn]   = np.array(t) - 100
            self.real['pegs'][stn]   = np.array(y)

class Specfem(PegsBase):
    def __init__(self, filepath, **kwargs):
        super().__init__(filepath=filepath, **kwargs)

        # Specfem needs to reverse Z acceleration relative to other codes.
        # So that PEGS = zacc + grad phi
        self.set_z_dir('zacc', -1)

        self.plot_colour = '#DDAA33'

class Qssp(PegsBase):
    def __init__(self, filepath, **kwargs):
        super().__init__(filepath=filepath, **kwargs)

        # Has down as positive:
        self.set_z_dir('grav', -1)
        self.set_z_dir('zacc', -1)

        self.plot_colour = '#004488'

class Axitra(PegsBase):
    def __init__(self, filepath, **kwargs):
        super().__init__(filepath=filepath, **kwargs)

        self.plot_colour = 'k'

    def _compute_pegs(self):
        # AXITRA data seems to be of different lengths
        for stn in self.stns:
            # Get shortest length:
            minlen = np.min([len(self.grav[stn]), len(self.zacc[stn])])
            self.time[stn]  = self.time[stn][:minlen]
            self.grav[stn]  = self.grav[stn][:minlen]
            self.zacc[stn]  = self.zacc[stn][:minlen]

        super()._compute_pegs()



def create_pegs(fname, code, scale=1e9):
    cu = code.upper()
    if   cu == 'SPECFEM':
        return Specfem(fname, scale=scale)
    elif cu == 'QSSP':
        return Qssp(fname, scale=scale)
    elif cu == 'AXITRA':
        return Axitra(fname, scale=scale)


