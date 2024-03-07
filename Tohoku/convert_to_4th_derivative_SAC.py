import matplotlib.pyplot as plt
import obspy
from obspy.core.stream import Stream
from obspy.core.trace import Trace
import numpy as np
from derivative.dglobal import Spectral
from scipy.signal.windows import tukey
from scipy.interpolate import CubicSpline


networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'XAN', 'ULN']

n = 10

dpath = "./forward_data_and_plots/elastic_test/elastic/"
fpath = f"../adjoint_source/4th_derivative/spline_deriv/n{n}"

rawdir = ''
outdir = 'Zacc'

for stn in stations:
    network = networks[stn]

    # Load data and get stats:
    if outdir=='Grav':
        d = np.loadtxt(f"{dpath}/FD_GRAV/{stn}_4.{network}.MXZC.PGRAV.sem.ascii")
    elif outdir=='Zacc':
        d = np.loadtxt(f"{dpath}/OUTPUT_FILES_ELASTIC/{stn}_0.{network}.MXZA.sem.ascii")
    else:
        raise ValueError('wrong value for outdir')

    time = d[:,0]
    acc  = d[:,1]

    dt = np.mean(d[1:, 0] - d[:-1, 0])


    # Conversion using FD:
    acc_2dx =  (-acc[4:] + 16*acc[3:-1] -30*acc[2:-2] + 16*acc[1:-3] - acc[:-4])/(12*dt*dt)

    # Conversion using spectral approach:
    #spec = Spectral()
    #acc_2dx = spec.compute_2nd(time, tukey(len(time), 0.15)*acc, None)[0,:]


    cs = CubicSpline(time[::n], acc[::n], bc_type='natural')

    #fig, ax = plt.subplots(2)
    #ax[0].plot(time, acc)
    #ax[0].plot(time[::n], acc[::n], 'x')
    #ax[0].plot(time, cs(time))
    #ax[1].plot(time, cs(time, 2))
    #ax[1].plot(time[2:-2], acc_2dx)
    #ax[1].axvline(198-105)
    #plt.show()
    acc_2dx = cs(time, 2)



    tr = Trace()
    tr.data = acc_2dx
    tr.stats.delta = dt



    st = Stream()
    st += tr
    st.write(filename=f'{fpath}/{outdir}/{stn}.sac')





