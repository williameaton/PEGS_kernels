import obspy
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np

networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'ULN'] # 'XAN',

station = 'MDJ'

fpath = f"../forward_simulations/Tohoku/data/globe_isotropic_prem/assuming_perfect_sphere/NEX256"

chl    = 'Z'


for stn in stations:
    network = networks[stn]

    # Load data and get stats:
    d = np.loadtxt(f"{fpath}/raw/{network}.{stn}.MX{chl}.sem.ascii")
    disp = d[:,1]

    dt = np.mean(d[1:, 0] - d[:-1, 0])

    # Compute second time derivative:
    acc = np.zeros(len(disp))

    acc[2:-2] =  (-disp[4:] + 16*disp[3:-1] - 30*disp[2:-2] + 16*disp[1:-3] - disp[:-4])/(12*dt*dt)
    acc[-2:] = acc[-3]

    tr = Trace()
    tr.data = acc
    tr.stats.delta = dt

    st = Stream()
    st += tr
    st.write(filename=f'{fpath}/{chl}acc/{stn}.sac')





