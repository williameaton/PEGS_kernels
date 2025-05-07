# Script to taper the forward wavefield and make it an adjoint source
import obspy
from obspy.core.trace import Trace
import numpy as np
import matplotlib.pyplot as plt

ttime= 20
starttaper = 320


# Load in the channels

stn = 'ULN_0'
net = 'IU'

fig, ax = plt.subplots(4)

ichlax = 1
for chl in ['Z', 'N', 'E']:


    # Load acceleration trace:
    trace = np.loadtxt(f'./OUTPUT_FILES/{stn}.{net}.MX{chl}V.sem.ascii')
    t = trace[:,0]
    d = trace[:,1]

    tr = Trace()

    tr.data = d
    tr.stats.delta = np.mean(t[1:] - t[:-1])
    tr = tr.filter('lowpass', freq=0.02)


    ax[ichlax].plot(t,d,'k')

    # Get filtered:
    d = tr.data
    ax[ichlax].plot(t,d,'b')


    # Create taper:

    taper = t*0 + 1
    tmax = t[-1]
    tmask =  np.logical_and(t>= starttaper, t<= starttaper + ttime)

    tupper = t[tmask]
    tmin = tupper[0]

    taper[tmask] = 0.5 + 0.5*np.cos( np.pi*(tupper - tmin)/(ttime) )

    taper[t>=starttaper+ttime] = 0

    tapered_trace = taper*d


    ax[0].plot(t, taper)

    ax[ichlax].plot(t, tapered_trace, 'r')


    ichlax += 1


    # Save the trace to SEM dir:
    np.savetxt(f'./SEM/{stn}.{net}.MX{chl}.adj', X=np.array([t,tapered_trace]).T )

plt.show()
