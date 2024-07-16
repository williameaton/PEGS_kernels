import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../classes')
from pegs import create_pegs


def obspy_gen_mpl(tr):
    x = np.linspace(0, tr.stats.npts*tr.stats.delta,  tr.stats.npts)
    y = tr.data
    return x,y

# Fiddle with order slightly so SPECFEM can plot on top
real_clr = '#BB5566'


# Some plotting stuff:
figsize = (6, 7.5)
lw_syn = 2
grid_lines_on = False
vlinmin = 0.015


# Load specfem and real data
spec = create_pegs('data/globe_isotropic_prem/assuming_perfect_sphere/NEX256', code='SPECFEM')
spec.load()
spec.load_real_data("./data/Vallee_data/real_data")
spec.cut_to_plot()

# Load QSSP
qssp = create_pegs('data/qssp/prem_elastic', code='QSSP')
qssp.load()
qssp.cut_to_plot()

# Load Axitra
axitra = create_pegs('data/Vallee_data/synthetic/', code='AXITRA', scale=1)
axitra.load()
axitra.cut_to_plot()



# Generate plots:
fig, ax = plt.subplots(figsize=figsize)
fig3D, ax3D = plt.subplots(figsize=figsize)
for f in [fig, fig3D]:
    f.set_tight_layout(True)

# Add grid lines
if grid_lines_on:
    for j in range(4):
        ax.axvline(j*100, ymin=vlinmin, ymax=1,  color='k',linestyle='--', alpha=0.05)


leg_lines = []

# Loop over each station
ictr = 0
for stn in spec.stns:

    # Plot real data from Vallee supplementary:
    for a in [ax, ax3D]:
        t = spec.real['time'][stn]
        p = spec.real['pegs'][stn]
        lreal, = a.plot(t[t < 0],  9 + p[t <  0] - ictr, color=real_clr, linewidth=1, alpha=0.2)
        lreal, = a.plot(t[t >= 0], 9 + p[t >= 0] - ictr, color=real_clr, linewidth=2, alpha=0.5)


    # Plot each code for this station
    for icode in [axitra, qssp, spec]:
        ltmp, = ax.plot(icode.time[stn], 9 + icode.pegs[stn] - ictr, color=icode.plot_colour, linewidth=lw_syn)

        # Only add to legend string once
        if ictr == 0:
            leg_lines.append(ltmp)

    # Add station name and horizontal lines
    for a in [ax, ax3D]:
        a.text(x=-147, y=9 - ictr - 0.1, s=stn, weight="bold")
        a.axhline(y=9 - ictr, xmin=0/450, xmax=100/450, color='grey', linewidth=1, alpha=0.3)

    ictr += 1


# Format legend
leg_lines.append(lreal)
legend = ax.legend(leg_lines, ['AXITRA','QSSP','SPECFEM', 'Real data'], fontsize=8)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)


for a in [ax, ax3D]:

    a.set_ylim([0, 9.5])
    a.set_xlim([-100, 350])
    a.axvline(0, ymin=vlinmin, ymax=1, linewidth=2, color='k')
    a.set_xlabel('Time since rupture initiation [s]', weight='bold')

    for xlabel_i in a.axes.get_yticklabels():
        xlabel_i.set_fontsize(10.0)
        xlabel_i.set_visible(False)


    a.errorbar(x=300, y=2.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
    a.text(x=305, y=2.4, s=r'1 $nm/s^2$' )


    a.tick_params(labelleft=False, left=False)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.spines['left'].set_visible(False)


    a.set_xticks([-100, 0, 100, 200, 300])
    a.set_xticklabels(('-100', '0', '100', '200', '300'))

#for path in ['./figures/', '../../Figures/pdfs/']:
#    fig.savefig(f'{path}/PEGS_tohoku_benchmark.pdf', format='pdf')
plt.show()

