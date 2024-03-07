import obspy
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np

networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'XAN', 'ULN']

fpath = "./adjoint_source/ULN_DATA/"

grav   =  False # false for acc
rawdir = 'OUTPUT_FILES'
chl    = 'E'

for grav in [True, False]:

    if grav:
        outdir = f'{chl}grav'
    else:
        outdir = f'{chl}acc'


    for stn in ['ULN']:
        network = networks[stn]

        # Load data and get stats:
        if grav:
            d = np.loadtxt(f"{fpath}/{rawdir}/{stn}_0.{network}.MX{chl}C.PGRAV.sem.ascii")
        else:
            d = np.loadtxt(f"{fpath}/{rawdir}/{stn}_0.{network}.MX{chl}A.sem.ascii")
        dt = np.mean(d[1:, 0] - d[:-1, 0])

        tr = Trace()
        tr.data = d[:,1]
        tr.stats.delta = dt

        st = Stream()
        st += tr
        st.write(filename=f'{fpath}/{outdir}/{stn}.sac')





