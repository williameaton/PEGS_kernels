import numpy as np
import matplotlib.pyplot as plt
import obspy
from colour_schemes import hex  # (pip install -i https://test.pypi.org/simple/ colour-schemes==0.0.1)

def obspy_gen_mpl(tr):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data
    return x,y

# Fiddle with order slightly so SPECFEM can plot on top
vibrant = hex.Hexes().hex['HighContrast_PT'][1::2][::-1]
real_clr ='#BB5566'


stn_list = ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO']
prem_cutoff = [57.95,   77.77,  163.81,  177.7 ,  237.43,  276.4 ,  292.3 , 343.65,  344. ]

fig, ax = plt.subplots(figsize=(6, 7.5))

leg_lines = []

# Some factors:
grid_lines_on = False
vlinmin = 0.015

lw_syn = 2

clrctr = 0
for src in [#'anelastic_prem_iso2c',
            'data/qssp/prem_elastic',
            'data/elastic_test/elastic', # 1D_isotropic_prem
            #'elastic_test/elastic_192',   # 1D_isotropic_prem
            #'S40RTS'
            ]:


    grav = []


    if grid_lines_on:
        for j in range(4):
            ax.axvline(j*100, ymin=vlinmin, ymax=1,  color='k',linestyle='--', alpha=0.05)

    ictr = 0
    for stn in stn_list:
        # Add station name
        ax.text(x=-135, y=9 - ictr - 0.1, s=stn, weight="bold")


        # Load the real data from Vallee Macro:
        if clrctr==0: #only do it once
            t, yG = obspy_gen_mpl(obspy.read(f'./data/Vallee_data/real_data/{stn}_acc_Z')[0])
            # Split linewidth around 1 = 0
            t = t - 100
            l0, =ax.plot(t[t<0], 9 + yG[t<0] - ictr, color=real_clr, linewidth=1, alpha=0.2)
            l0, =ax.plot(t[t>=0], 9 + yG[t>=0] - ictr, color=real_clr, linewidth=2, alpha=0.5)

        # Plot the Vallee synthetics:
        if clrctr==0: #only do it once
            tG,yG = obspy_gen_mpl(obspy.read(f'./data/Vallee_data/SAC_deltaG/{stn}_Z.SAC')[0])
            t,y   = obspy_gen_mpl(obspy.read(f'./data/Vallee_data/SAC_AccInduit/{stn}_Z.sac')[0])
            m = np.min([len(yG), len(y)])
            lv, = ax.plot(t[:m], 9 + yG[:m] - y[:m]  - ictr, 'k', linewidth=lw_syn)


        # Get time cutoff
        cut = prem_cutoff[8-ictr]

        # Add background lines
        if grid_lines_on:
            ax.axhline(ictr, color='k', linestyle='--', alpha=0.05)

        # Load gravity data:
        t,yG = obspy_gen_mpl(obspy.read(f'{src}/Grav/{stn}_proc.sac')[0])
        yG*=-1e9

        m = t < cut + 35
        if 'qssp' in src:
            t = t[m]
            yG = -yG[m]
        else:
            t = t[m] - 35
            yG = yG[m]
        #lgrav, = ax.plot(t, 9 + yG - ictr, 'teal', linewidth=1.3, alpha=0.5)




        # Plot the acceleration data
        t,y = obspy_gen_mpl(obspy.read(f'{src}/Zacc/{stn}_proc.sac')[0])
        y*=-1e9

        m = t < cut + 35

        if 'qssp' in src:
            t = t[m]
        else:
            t = t[m] - 35

        y = y[m]
        #lacc, = ax.plot(t, 9 - y - ictr, 'r', linewidth=1.3, alpha=0.5)



        # Get shortest len:
        minlen = np.min([len(y), len(yG)])
        t  = t[:minlen]
        yG = yG[:minlen]
        y = y[:minlen]

        # PLOT PEGS!
        l3, = ax.plot(t[t>0], 9 -yG[t>0] - y[t>0]  - ictr, vibrant[clrctr], linewidth=lw_syn)
        l4 = ax.axhline(y=9 - ictr, xmin=0/450, xmax=100/450, color=vibrant[clrctr], linewidth=1, alpha=0.3)



        pegs = np.array([t,  -yG - y]).T
        #np.savetxt(fname=f'{src}/PEGS/{stn}_pegs', X=pegs)


        ax.set_ylim([0, 9.5])
        ictr += 1
    # append the synthetic
    leg_lines.append(l3)
    clrctr += 1

# Add in the vallee synthetics and Real line
leg_lines.append(lv)
leg_lines.append(l0)


#legend = ax.legend(leg_lines, ['Vallée, et al., 2017 data', 'PREM Mw 8.5', 'PREM Mw 9.0', 'Real acc.'])
legend = ax.legend(leg_lines, ['QSSP','SPECFEM', 'Vallée, et al., 2017', 'Real acc.'], fontsize=8)
#legend = ax.legend(leg_lines, ['PREM', 'Vallée, et al., 2017', 'Real acc.'])
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)


ax.set_xlim([-100,350])
ax.set_xlabel('Time [s]', weight='bold')


ax.axvline(0, ymin=vlinmin, ymax=1, linewidth=2, color='k')

for xlabel_i in ax.axes.get_yticklabels():
    xlabel_i.set_fontsize(10.0)
    xlabel_i.set_visible(False)


ax.errorbar(x=300, y=2.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
ax.text(x=305, y=2.4, s=r'1 $nm/s^2$' )


ax.tick_params(labelleft=False, left=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)


ax.set_xticks([-100, 0, 100, 200, 300])
ax.set_xticklabels(('-100', '0', '100', '200', '300'))

plt.show()