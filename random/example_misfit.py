import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss


import numpy as np
import matplotlib.pyplot as plt
import obspy
from wetools import obspy_gen_mpl


clr_key = {'real': '#BB5566'}

clrs = ['#DDAA33' ,'#004488']

stn = 'INU'

fig, ax = plt.subplots(figsize=(8, 5))

leg_lines = []




M = 1000
w = 3.5

w1 = ss.morlet(M, 1, w).imag
w2 = ss.morlet(M, 2, w).imag
t  = np.linspace(0, 80, M)

ax.plot(t,w1, clrs[0], linewidth=3)
ax.plot(t,w2,  clrs[1], linewidth=3)



ax.set_xlim([-2, 80])
ax.set_ylim([-0.7, 0.7])

ax.set_xticks([0, 20, 40, 60, 80])
ax.set_xticklabels(('0', '20', '40', '60', '80'))

ax.set_yticks([-0.5, 0, 0.5])
ax.set_yticklabels(('-1', '0', '1'))

fs = 12
ax.set_xlabel('Time [s]', weight='bold', fontsize=fs)
ax.set_ylabel('Vertical acceleration [nm/s/s]', weight='bold', fontsize=fs)


#ax.tick_params(labelleft=False, left=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)




legend = ax.legend([r'Synthetic', 'Real'], fontsize=13)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)


ax.text(s='STATION: INU', x=0, y=-1, fontsize=11, weight='bold')

plt.show()





