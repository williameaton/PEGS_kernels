import numpy as np
import matplotlib.pyplot as plt
import obspy
from wetools import obspy_gen_mpl


stn_list = ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO']
prem_cutoff = [57.95,   77.77,  163.81,  177.7 ,  237.43,  276.4 ,  292.3 , 343.65,  344. ]

fig, ax = plt.subplots(figsize=(5, 7.5))

leg_lines = []

clrs = ['teal', 'orange', 'red', 'olive']
clrctr = 0
for src in ['anelastic_prem_iso2c',
            #'elastic_test/anelastic', # 1D_isotropic_prem
            'elastic_test/elastic_192',   # 1D_isotropic_prem
            'S40RTS']:

    datadirs = [ f'./specfemx/{src}']


    grav = []

    ax.set_title('SPECFEMX Tohoku 2011 PEGS signals')

    for j in range(4):
        ax.axvline(j*100,  color='k',linestyle='--', alpha=0.05)

    ictr = 0
    for stn in stn_list:

        # Load the real data from Vallee Macro:
        if clrctr==0: #only do it once
            t, yG = obspy_gen_mpl(obspy.read(f'./Vallee_data/real_data/{stn}_acc_Z')[0])
            # Split linewidth around 1 = 0
            t = t - 100
            l0, =ax.plot(t[t<0], 9 + yG[t<0] - ictr, 'r', linewidth=1, alpha=0.5)
            l0, =ax.plot(t[t>=0], 9 + yG[t>=0] - ictr, 'r', linewidth=2, alpha=0.5)

        # Plot the Vallee synthetics:
        if clrctr==0: #only do it once
            tG,yG = obspy_gen_mpl(obspy.read(f'./Vallee_data/SAC_deltaG/{stn}_Z.SAC')[0])
            t,y   = obspy_gen_mpl(obspy.read(f'./Vallee_data/SAC_AccInduit/{stn}_Z.sac')[0])
            m = np.min([len(yG), len(y)])
            lv, = ax.plot(t[:m], 9 + yG[:m] - y[:m]  - ictr, 'k-', linewidth=1.5)

        cut = prem_cutoff[8-ictr]
        # Add background lines
        ax.axhline(ictr, color='k',linestyle='--', alpha=0.05)

        # Load gravity data:
        t,yG = obspy_gen_mpl(obspy.read(f'{src}/Grav/{stn}_proc.sac')[0])
        yG*=-1e9

        m = t < cut + 35
        t = t[m] - 35
        yG = yG[m]
        #l1, =ax.plot(t, 9 + yG - ictr, 'g', linewidth=1.3, alpha=0.5)

        # Add station name
        ax.text(x=-150, y=9 - ictr - 0.1, s=stn)

        # Plot the acceleration data
        t,y = obspy_gen_mpl(obspy.read(f'{src}/Zacc/{stn}_proc.sac')[0])
        y*=-1e9

        m = t <  cut + 35
        t = t[m] - 35
        y = y[m]
        #l2, = ax.plot(t, 9 - y - ictr, 'b', linewidth=1.3, alpha=0.5)


        # Get shortest len:
        minlen = np.min([len(y), len(yG)])
        t  = t[:minlen]
        yG = yG[:minlen]
        y = y[:minlen]

        # PLOT PEGS!
        #l3, = ax.plot(t, 9  - y - ictr, clrs[clrctr], linewidth=1)
        #l3, = ax.plot(t, 9  - yG - ictr, clrs[clrctr], linewidth=1)
        l3, = ax.plot(t[t<0], 9 -yG[t<0] - y[t<0]  - ictr, clrs[clrctr], linewidth=1)
        l3, = ax.plot(t[t>=0], 9 -yG[t>=0] - y[t>=0]  - ictr, clrs[clrctr], linewidth=1.5)
        l4 = ax.axhline(y=9 - ictr, xmin=0/450, xmax=100/450, color=clrs[clrctr])


        pegs = np.array([t,  -yG - y]).T
        np.savetxt(fname=f'{src}/PEGS/{stn}_pegs', X=pegs)


        ax.set_ylim([0, 9.5])
        ictr += 1
    # append the synthetic
    leg_lines.append(l3)
    clrctr += 1

# Add in the vallee synthetics and Real line
leg_lines.append(lv)
leg_lines.append(l0)

legend = ax.legend(leg_lines, ['PREM', 'S40RTS', 'Vallée, et al., 2017', 'Real acc.'])
#legend = ax.legend(leg_lines, ['PREM', 'Vallée, et al., 2017', 'Real acc.'])
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)


ax.set_xlim([-100,350])
ax.set_xlabel('Time [s]')

ax.spines[['right', 'top']].set_visible(False)

ax.axvline(0, linewidth=2, color='k')

for xlabel_i in ax.axes.get_yticklabels():
    xlabel_i.set_fontsize(0.0)
    xlabel_i.set_visible(False)
for tick in ax.axes.get_yticklines():
    tick.set_visible(False)

ax.errorbar(x=300, y=2.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
ax.text(x=305, y=2.4, s=r'1 $nm/s^2$' )

plt.show()