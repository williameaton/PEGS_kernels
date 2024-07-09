import numpy as np
import matplotlib.pyplot as plt
import obspy
import matplotlib as mpl
from matplotlib.patches import Rectangle
from wetools import obspy_gen_mpl

clr_key = {'real': '#BB5566'}

clrs = ['#DDAA33' ,'#004488']

stn = 'INU'

fig, ax = plt.subplots(figsize=(8, 5))

leg_lines = []



src = 'anelastic_prem_iso2c'


# Load the real data from Vallee Macro:
t, yG = obspy_gen_mpl(obspy.read(f'./data/Vallee_data/real_data/{stn}_acc_Z')[0])
# Split linewidth around 1 = 0
t = t - 100
l0, =ax.plot(t[t>=0], yG[t>=0], color=clr_key['real'], linewidth=3)



cut = 77.77


# Load gravity data:
t,yG = obspy_gen_mpl(obspy.read(f'./data/{src}/Grav/{stn}_proc.sac')[0])
yG*=-1e9

m = t < cut + 35
t = t[m] - 35
yG = yG[m]



# Plot the acceleration data
t,y = obspy_gen_mpl(obspy.read(f'./data/{src}/Zacc/{stn}_proc.sac')[0])
y*=-1e9

m = t <  cut + 35
t = t[m] - 35
y = y[m]


# Get shortest len:
minlen = np.min([len(y), len(yG)])
t  = t[:minlen]
yG = yG[:minlen]
y = y[:minlen]



# PLOT PEGS!
l3, = ax.plot(t[t>0],  y[t>0] , clrs[0], linewidth=3)
l3, = ax.plot(t[t>0],  - yG[t>0] , clrs[1], linewidth=3)
l3, = ax.plot(t[t>0], -yG[t>0] - y[t>0]  , 'k', linewidth=3)



ax.set_xlim([-2, 80])
ax.set_ylim([-1.2, 3.2])

ax.set_xticks([0, 20, 40, 60, 80])
ax.set_xticklabels(('0', '20', '40', '60', '80'))

ax.set_yticks([-1, 0, 1, 2, 3])
ax.set_yticklabels(('-1', '0', '1', '2', '3'))

fs = 12
ax.set_xlabel('Time [s]', weight='bold', fontsize=fs)
ax.set_ylabel('Vertical acceleration [nm/s/s]', weight='bold', fontsize=fs)


#ax.tick_params(labelleft=False, left=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)




legend = ax.legend([r'Real acc. ($A$)', 'Ground acc. ($\partial_t^2 \mathbf{s}$)', 'Gravity perturb. ($\mathbf{\delta g}$)', 'Net PEGS signal ($\mathbf{\delta g} - \partial_t^2 \mathbf{s}$)'], fontsize=13)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)


ax.text(s='STATION: INU', x=0, y=-1, fontsize=11, weight='bold')
ax.text(s='M9.0 TŌHOKU-OKI, JAPAN (2011) ', x=40, y=-1, fontsize=11, weight='bold')

plt.show()