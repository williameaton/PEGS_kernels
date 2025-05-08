import numpy as np
import matplotlib.pyplot as plt
from copy import copy
networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }

mdir = 'forward_data_and_plots/elastic_test/elastic_192'
stn  = 'INCN'
network = networks[stn]

# Load the PEGS signal:
pgs = np.loadtxt(f"{mdir}/PEGS/{stn}_pegs")
tPGS = pgs[:,0] - 70
dPGS = pgs[:,1]


# Load a trace to get the timings necessary:
tr = np.loadtxt(f"{mdir}/OUTPUT_FILES/{stn}_0.{network}.MXZA.sem.ascii")
ttr = tr[:,0]



# Resample the PEGS signal at correct rate:
old_dt = np.mean(tPGS[1:] - tPGS[:-1])

new_dt = np.mean(tr[1:,0] - tr[:-1,0])

# Create resampled data to fill in:
new_d = copy(ttr)*0

try:
    for i in range(1,len(new_d)):
        t_i = ttr[i]

        gap = tPGS[i] - tPGS[i-1]

        val = dPGS[i-1] + ((dPGS[i] - dPGS[i-1])/gap) *  (t_i - tPGS[i-1])

        new_d[i] = val
except:
    pass

new_d[len(dPGS)] = dPGS[-1]




fig, ax = plt.subplots(figsize=(12,7.5))

#ax.plot(tPGS, dPGS)

# PLOT RESAMPLED PEGS
#ax.plot(ttr, new_d, 'k')

# PLOT RESAMPLED PEGS UP TO POINT WHERE P WAVE ARRIVES
ax.plot(ttr[ttr<107.7], new_d[ttr<107.7], 'k')



# Now we need to taper the bloody thing!


ttime= 2.5
starttaper = tPGS[-1] - ttime

# Create taper:

taper = ttr * 0 + 1
tmax = ttr[-1]
tmask = np.logical_and(ttr >= starttaper, ttr <= starttaper + ttime)

tupper = ttr[tmask]
tmin = tupper[0]

taper[tmask] = 0.5 + 0.5 * np.cos(np.pi * (tupper - tmin) / (ttime))

taper[ttr >= starttaper + ttime] = 0

tapered_trace = taper * new_d


#ax.plot(ttr, tapered_trace, 'r')

ax.set_xlabel('Time [s]')
ax.set_ylabel('PEGS Acceleration [nm/s^2]')


np.savetxt(f'{mdir}/ADJ/{stn}.{network}.MXZ.adj', X=np.array([ttr, tapered_trace]).T)

ax.legend(['Original PEGS acceleration', 'Tapered trace'])

ax.set_xlim([-105, 194])

plt.show()