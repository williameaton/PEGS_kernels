import numpy as np
import matplotlib.pyplot as plt
import obspy
from wetools import obspy_gen_mpl


prem_cutoff = {'MAJO':57.95,   'INU':77.77,  'MDJ':163.81,  'INCN':177.7 ,
               'NE93':237.43,  'BJT':276.4 , 'MA2':292.3 ,  'XAN':343.65, 'ULN': 344. }

stn = 'INU'
cut = prem_cutoff[stn]

fig, ax = plt.subplots(3, sharex=True, figsize=(9.5,6.5))

# Real data:
tR, yR = obspy_gen_mpl(obspy.read(f'./forward_data_and_plots/Vallee_data/real_data/{stn}_acc_Z')[0])
tR = tR - 100
ax[0].plot(tR[tR < 0],yR[tR<0], alpha=0.5, color='r')
m = tR>=0
tR = tR[m]
yR = yR[m]
ax[0].plot(tR,yR, 'r', alpha=1)


pegs = np.loadtxt(f"./forward_data_and_plots/anelastic_prem_iso2c/PEGS/{stn}_pegs")
m =  pegs[:,0]>=0
tP = pegs[m,0]
yP = pegs[m,1]
ax[0].plot(tP,yP, 'k')


# Reduce to shortest len:
minlen = np.min([len(tP), len(tR)])
yP = yP[:minlen]
tP = tP[:minlen]
yR = yR[:minlen]
tR = tR[:minlen]


# LOG MISFIT:
log_chi =0.5 *  (np.log(np.abs(yP)/ np.abs(yR))  )**2

# Plot individual logs
l0, = ax[1].plot(tP, np.log(np.abs(yP)), 'k')
l1, = ax[1].plot(tR, np.log(np.abs(yR)), 'r')
ax[1].set_ylabel('Log acceleration')

axchi = ax[1].twinx()
l2, = axchi.plot(tR, log_chi, 'b')
axchi.set_ylabel('Log Misfit')

# L2 MISFIT:
ax[2].plot(tP,  0.5 * ((yP - yR)**2), 'b')
ax[2].set_ylabel('L2 misfit')

ax[1].legend([l0,l1,l2],['Synthetic PEGS', 'Real data', 'Log misfit'], loc='upper left')

ax[-1].set_xlabel( 'Time after rupture [s]')
ax[-1].set_xlim([-100, cut])

fig.suptitle(f'Station: {stn}')

for i in range(3):
    ax[i].axvline(0, color='k', linewidth=1.5)


plt.savefig(fname=f"{stn}_misfit.pdf", format='pdf')