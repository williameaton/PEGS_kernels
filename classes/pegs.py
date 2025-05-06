from abc import ABC, abstractmethod
from wetools.funcs import obspy_gen_mpl
import numpy as np
import obspy

class PegsBase(ABC):

    @abstractmethod
    def __init__(self, filepath, **kwargs):

        # Default kwargs:
        defaultKwargs = {'stations'   : ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO'],
                         'cutoff'     : {'MAJO': 57.95,  'INU':  77.77, 'MDJ':  163.81,  'INCN': 177.7 ,  'NE93': 237.43, 'BJT':  276.4 , 'MA2':  292.3 ,'XAN':  343.65,  'ULN': 344. },
                         'scale'      : 1e9,
                         'three_component': False,}
        kwargs = defaultKwargs | kwargs

        # Store inputs
        self.fpath = filepath

        # Possible kwargs:
        self.stns  = kwargs['stations']
        self.scale = kwargs['scale']
        self.cutoff = kwargs['cutoff']
        self.three_component = kwargs['three_component']

        # Setup
        self.time = {}
        self.real = {'time' : {},
                     'pegs' : {}}

        self.z_pegs = {}
        self.z_grav = {}
        self.z_acc  = {}

        self.traces = {'time'  : self.time,
                       'z_pegs': self.z_pegs,
                       'z_grav': self.z_grav,
                       'z_acc' : self.z_acc}

        self.ncomponents = 1
        self.channels = 'z'

        # Add the extra components
        if self.three_component:
            self.n_pegs = {}
            self.n_grav = {}
            self.n_acc  = {}

            self.e_pegs = {}
            self.e_grav = {}
            self.e_acc  = {}
            self.traces['e_pegs'] =  self.e_pegs
            self.traces['e_grav'] =  self.e_grav
            self.traces['e_acc']  =  self.e_acc

            self.traces['n_pegs'] =  self.n_pegs
            self.traces['n_grav'] =  self.n_grav
            self.traces['n_acc']  =  self.n_acc

            self.ncomponents = 3
            self.channels    = 'zne'


        self.z_dir = {'grav': 1,
                      'acc': 1}

        self.nstns = len(self.stns)




    def load(self):

        # Load z data
        self._load_data('grav')
        self._load_data('acc')

        if self.three_component:
            for chl in 'en':
                self._load_data('grav', chl=chl)
                self._load_data('acc', chl=chl)

        self._compute_pegs()

    def _load_data(self, dtype, chl='z'):

        supchl = chl.upper()

        for stn in self.stns:
            # Load to numpy arrays
            t, y = obspy_gen_mpl(obspy.read(f'{self.fpath}/{supchl}{dtype}/{stn}_proc.sac')[0])

            # Add to relevant array:
            self.traces['time'][stn] = np.array(t)

            # Deals with direction for positive
            if chl=='z':
                mult = self.z_dir[dtype]
            else:
                mult = 1

            self.traces[f"{chl}_{dtype}"][stn]  = np.array(y*mult) * self.scale


    def _compute_pegs(self):
        for chl in self.channels:
            for stn in self.stns:
                self.traces[f"{chl}_pegs"][stn] = self.traces[f"{chl}_grav"][stn] - self.traces[f"{chl}_acc"][stn]


    def set_z_dir(self, dtype, value):
        if np.abs(value)!=1:
            raise ValueError('value must be +-1')
        if dtype!="grav" and dtype!="acc":
            raise ValueError(f"dtype must be 'grav' or 'acc'")

        self.z_dir[dtype] = value

    def cut_to_plot(self, plusval=35):
        for chl in self.channels:
            for stn in self.stns:
                t = self.time[stn]
                m = t < self.cutoff[stn] + plusval

                self.time[stn] = t[m]

                self.traces[f"{chl}_grav"][stn] = self.traces[f"{chl}_grav"][stn][m]
                self.traces[f"{chl}_acc" ][stn] = self.traces[f"{chl}_acc" ][stn][m]


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
        self.set_z_dir('acc', -1)

        self.plot_colour = '#DDAA33'
        self.plot_label  = 'SPECFEM'

class Qssp(PegsBase):
    def __init__(self, filepath, **kwargs):
        super().__init__(filepath=filepath, **kwargs)

        # Has down as positive:
        self.set_z_dir('grav', -1)
        self.set_z_dir('acc',  -1)

        self.plot_colour = '#004488'
        self.plot_label  = 'QSSP'

class Axitra(PegsBase):
    def __init__(self, filepath, **kwargs):
        super().__init__(filepath=filepath, **kwargs)

        self.plot_colour = '#BB5566'
        self.plot_label  = 'AXITRA'
    def _compute_pegs(self):
        # AXITRA data seems to be of different lengths
        for stn in self.stns:
            # Get shortest length:
            minlen = np.min([len(self.z_grav[stn]), len(self.z_acc[stn])])
            self.time[stn]   = self.time[stn][:minlen]
            self.z_grav[stn]  = self.z_grav[stn][:minlen]
            self.z_acc[stn]   = self.z_acc[stn][:minlen]

        super()._compute_pegs()



def create_pegs(fname, code, scale=1e9, three_component=False):
    cu = code.upper()
    if   cu == 'SPECFEM':
        return Specfem(fname, scale=scale, three_component=three_component)
    elif cu == 'QSSP':
        return Qssp(fname, scale=scale, three_component=three_component)
    elif cu == 'AXITRA':
        return Axitra(fname, scale=scale, three_component=three_component)


