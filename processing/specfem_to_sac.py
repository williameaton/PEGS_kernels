# Converts SPECFEM ascii output to SAC
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np

# SPECIFY NETWORKS AND STATIONS:
# todo: modify to read specfem STATIONS file
networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'XAN', 'ULN']

# Will load ascii from ```path/rawdir``` and save sac files in ```path/{chl}grav and``` ```path/{chl}acc```
fpath = f"/Users/eaton/Documents/Princeton/PEGS_kernels/data/MDJ_pegs_kernel_Zonly/forward_output/"


rawdir        = 'raw'   # SPECFEM ascii data directory
chl           = 'Z'     # Channel of interet
SPECFEM_GLOBE = False   # Set to true for spfm_globe ascii formatting

for grav in [True, False]:
    if grav:
        outdir = f'{chl}grav'
    else:
        outdir = f'{chl}acc'


    for stn in stations:
        network = networks[stn]

        # Load data and get stats:
        if grav:
            if SPECFEM_GLOBE:
                d = np.loadtxt(f"{fpath}/raw/{network}.{stn}.MX{chl}.C.PGRAV.sem.ascii")
            else:
                d = np.loadtxt(f"{fpath}/{rawdir}/{stn}.{network}.MX{chl}C.PGRAV.sem.ascii")
        else:
            if SPECFEM_GLOBE:
                d = np.loadtxt(f"{fpath}/raw/{network}.{stn}.MX{chl}A.sem.ascii")
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





