import matplotlib.pyplot as plt
import obspy
from scipy.signal.windows import tukey
from derivative.dglobal import Spectral
from wetools import obspy_gen_mpl
import numpy as np

# Load the pegs data:
stn = 'MDJ'


zacc_tr = obspy.read(f'../Tohoku/forward_data_and_plots/elastic_test/anelastic/not_cutoff/Zacc/{stn}_proc.sac')[0]
grav_tr = obspy.read(f'../Tohoku/forward_data_and_plots/elastic_test/anelastic/not_cutoff/Grav/{stn}_proc.sac')[0]

pegs_tr = zacc_tr.copy()
pegs_tr.data = zacc_tr.data + grav_tr.data

ntime = 1000
# Convert to numpy
time, zacc = obspy_gen_mpl(zacc_tr)
time, grav = obspy_gen_mpl(grav_tr)
time, pegs = obspy_gen_mpl(pegs_tr)


tmask = time < ntime
time = time[tmask]
zacc = zacc[tmask]
grav = grav[tmask]
pegs = pegs[tmask]

window = tukey(len(time), 0.1)
pegs *= window

fig, ax = plt.subplots(4)

ax[0].plot(time, zacc)
ax[1].plot(time, grav)

PEGS = zacc + grav

ax[2].plot(time, PEGS)
ax[2].plot(time, pegs)



# Spectal acceleration:
spec = Spectral()
ff = spec.compute_2nd(time,pegs, None)
ax[3].plot(time, ff[0,:])

# FD acceleration:
dx = np.mean(time[1:] - time[:-1])
acc_2dx = (-pegs[4:] + 16 * pegs[3:-1] - 30 * pegs[2:-2] + 16 * pegs[1:-3] - pegs[:-4]) / (12 * dx * dx)

ax[-1].plot(time[2:-2], acc_2dx)
plt.show()

