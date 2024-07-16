import numpy as np
import matplotlib.pyplot as plt
import obspy
import matplotlib as mpl
from matplotlib.patches import Rectangle
from wetools import obspy_gen_mpl
import sys
sys.path.append('../../classes')
from pegs import create_pegs


clr_key = {'real': '#BB5566'}
clrs    = ['#DDAA33' ,'#004488']
stn       = 'MDJ'


fig, ax = plt.subplots(figsize=(8, 5))

leg_lines = []

spec = create_pegs('data/globe_isotropic_prem/assuming_perfect_sphere/NEX256', code='SPECFEM')
spec.load()
spec.load_real_data("./data/Vallee_data/real_data")
spec.cut_to_plot(plusval=35)


# plot real data
l0, = ax.plot(spec.real['time'][stn], spec.real['pegs'][stn], color=clr_key['real'], linewidth=3)

# plot synthetics
l3, = ax.plot(spec.time[stn],  spec.grav[stn], clrs[0], linewidth=3)
l3, = ax.plot(spec.time[stn],  spec.zacc[stn], clrs[1], linewidth=3)
l3, = ax.plot(spec.time[stn],  spec.pegs[stn], 'k',     linewidth=3)



ax.set_xlim([-20, 180])
ax.set_ylim([-1.75, 1])

#ax.set_xticks([0, 20, 40, 60, 80])
#ax.set_xticklabels(('0', '20', '40', '60', '80'))

ax.set_yticks([-1.5, -1, -0.5, 0, 0.5])
ax.set_yticklabels(('-1.5', '-1', '-0.5', '0', '0.5'))

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
ax.text(s='M9.0 TŌHOKU-OKI, JAPAN (2011) ', x=-17, y=0.9, fontsize=11, weight='bold')

plt.show()