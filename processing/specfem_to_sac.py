import obspy
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np

networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'XAN', 'ULN']

station = 'MDJ'

fpath = f"../kernels/adjoint_source/SEPT_5"

rawdir = 'OUTPUT_FILES'
chl    = 'Z'

for grav in [True, False]:

    if grav:
        outdir = f'{chl}grav'
    else:
        outdir = f'{chl}acc'


    for stn in stations: #[f'{station}']:
        network = networks[stn]

        # Load data and get stats:
        if grav:
            d = np.loadtxt(f"{fpath}/{rawdir}/{stn}.{network}.MX{chl}C.PGRAV.sem.ascii")
        else:
            d = np.loadtxt(f"{fpath}/{rawdir}/{stn}.{network}.MX{chl}A.sem.ascii")
        dt = np.mean(d[1:, 0] - d[:-1, 0])
        print(dt)
        tr = Trace()
        tr.data = d[:,1]
        tr.stats.delta = dt

        st = Stream()
        st += tr
        outfname = f'{fpath}/{outdir}/{stn}.sac'
        st.write(filename=outfname)
        print(f"Written to {outfname}")





