import obspy
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np

networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2','XAN', 'ULN']

station = 'MDJ'

fpath = f"../forward_simulations/Tohoku/data/globe_isotropic_prem/no_perfect_sphere/NEX256"

chl    = 'Z'

for grav in [False, True]:

    if grav:
        outdir = f'Grav'
    else:
        outdir = f'{chl}acc'


    for stn in stations: #[f'{station}']:
        network = networks[stn]

        # Load data and get stats:
        if grav:
            d = np.loadtxt(f"{fpath}/raw/{network}.{stn}.MX{chl}.C.PGRAV.sem.ascii")
        else:
            d = np.loadtxt(f"{fpath}/raw/{network}.{stn}.MX{chl}A.sem.ascii")
        dt = np.mean(d[1:, 0] - d[:-1, 0])
        print(dt)
        tr = Trace()
        tr.data = d[:,1]
        tr.stats.delta = dt

        st = Stream()
        st += tr
        st.write(filename=f'{fpath}/{outdir}/{stn}.sac')





