# Plots Tohoku PEGS from SPECFEM compared with QSSP, AXITRA and the real data
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../classes')
from pegs import create_pegs
from scipy import signal
import pickle

# Some plotting stuff:
figsize_single  = (6, 7)         # Figure size
figsize_all  = (13, 7)         # Figure size
lw_syn = 1.5                  # Synthetic data line width
# grid_lines_on = True       # Grid lines?
vlinmin = 0.015             # Min y value of vertical lines
network = {'MAJO': 'IU', 'FUK': 'BO', 'SHR': 'BO', 'INU': 'G', 'BJT': 'IC', 'MDJ': 'IC', 'XAN': 'IC', 'INCN': 'IU',
           'MA2': 'IU', 'ULN': 'IU', 'NE93': 'YP'}

alpha3Dcodes = 1          # Alpha for codes other than 3D specfem in 3D plot

# Test plot options:
PT_bright   = ['#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377',]
PT_highcont_wwine = ['#DDAA33', '#BB5566', '#004488', '#882255']
PT_highcont_wteal = ['#DDAA33', '#BB5566', '#004488', '#009988']
PT_vibrant    = ['#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377']
PT_medcont    = ['#EECC66', '#EE99AA', '#6699CC', '#997700', '#994455', '#004488']
carto_antique = ["#855C75", "#D9AF6B", "#AF6458", "#736F4C", "#526A83", "#625377", "#68855C", "#9C9C5E", "#A06177", "#8C785D", "#467378", "#7C7C7C"]
carto_bold    = ["#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695", "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99"]
carto_prism   = ["#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05", "#CC503E", "#94346E", "#6F4070", "#994E95", "#666666"]
carto_safe    = ["#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"]
carto_vivid   = ["#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99"]

cbar = carto_antique


fig, ax = plt.subplots()
for i in range(12):
    ax.axhline(i, color=cbar[i], linewidth=5)
#plt.show()

# Original
#clrs = [cbar[0], cbar[1], cbar[2], cbar[3]]
#real_clr = "k"

# Safe version
#clrs = [cbar[1], cbar[8], cbar[2], cbar[4]]
#real_clr = cbar[11]

# Vivid version
#clrs = [cbar[5], cbar[8], cbar[1], cbar[6]]
#real_clr = cbar[11]


data_dir = '../../data/data_forward_pegs'


# Load the epicentral distances in km for annotation:
with open('../../Figures/data_for_plots/epi_dist.pkl', 'rb') as fp:
    epidists = pickle.load(fp)

# Load specfem for perfect sphere
spec = create_pegs(f'{data_dir}/globe_isotropic_prem/assuming_perfect_sphere/NEX256/', code='SPECFEM')
spec.load()

# 3D SPECFEM data:
spec3D = create_pegs(f'{data_dir}/S40RTS/', code='SPECFEM')
spec3D.load()
for stn in spec3D.stns:
    spec3D.time[stn] -= 35

spec3D.plot_label = 'S40RTS'

# Load real data
spec.load_real_data(f"./{data_dir}/Vallee_data/real_data")
spec.cut_to_plot()

# Load QSSP
qssp = create_pegs(f'{data_dir}/qssp/original_AK135/output/', code='QSSP')
qssp.load()
qssp.cut_to_plot()

# Load Axitra
axitra = create_pegs(f'{data_dir}/Vallee_data/synthetic/', code='AXITRA', scale=1)
axitra.load()
axitra.cut_to_plot()


# Antique version
#clrs = [cbar[0], , cbar[5], cbar[1], cbar[2]]

real_clr = carto_antique[9]
axitra.plot_colour = "#004488"
mineos_plot_colour = carto_antique[10]
qssp.plot_colour   = "#BB5566"
spec.plot_colour   = "#DDAA44"
spec3D.plot_colour = 'k'



# Create figures and set tight layouts:
fig, ax     = plt.subplots(1,3, figsize=figsize_all)     # Figure for 1D Prem
fig3D, ax3D = plt.subplots(figsize=figsize_single)     # Figure for 3D mantle simulation
all_axes = np.concatenate(([ax3D], ax))

for f in [fig, fig3D]:
    f.set_tight_layout(True)

# Add grid lines
# if grid_lines_on:
#     for j in range(4):
#         aax = np.concatenate(ax, ax3D)
#         for a in aax:
#             a.axvline(j*100, ymin=vlinmin, ymax=1,  color='k',linestyle='--', alpha=0.05)

# Stores line objects for mpl legend in loop
leg_lines = []
leg_lines_3D = []

# Loop over each station
ictr = 0
for stn in spec.stns:

    # Load real data from Vallee supplementary:
    t = spec.real['time'][stn]
    p = spec.real['pegs'][stn]

    # Using two lines for the 'before' and 'after' t=0 since they have different line thickness and transparency
    lreal, = ax[-1].plot(t[t < 0],  9 + p[t <  0] - ictr, color=real_clr, linewidth=1, alpha=0.2, label='Real Data')
    lreal, = ax[-1].plot(t[t >= 0], 9 + p[t >= 0] - ictr, color=real_clr, linewidth=2, alpha=0.2, label='Real Data')

    # Change transparency for 3D plot
    lreal, = ax3D.plot(t[t < 0],  9 + p[t <  0] - ictr, color=real_clr, linewidth=1, alpha=0.2, label='Real Data')
    lreal, = ax3D.plot(t[t >= 0], 9 + p[t >= 0] - ictr, color=real_clr, linewidth=2, alpha=0.2, label='Real Data')


    # Decided to add the Normal Modes Juhel code data:
    # This code was provided by Kevin Juhel
    # set sampling rate
    # set PEGS and time vector
    nm_pegs = np.load(f'{data_dir}/Juhel_2019_NM/KJ19/PEGS.' + network[stn] + '.' + stn + '.Z.pegs_proc.npy')
    time = np.arange(0, nm_pegs.size, 1.0) - 100
    timemask = np.logical_and(time < spec.cutoff[stn], time >=0)

    # Avoid a loop here over [ax ax3D] for duplicate label issues in legend
    lnm,   = ax[-1].plot(time[timemask], 9 - ictr + (10 ** 9 * nm_pegs[timemask]), color=mineos_plot_colour, linewidth=lw_syn, label='MINEOS')
    lnm3d, = ax3D.plot(time[timemask], 9 - ictr + (10 ** 9 * nm_pegs[timemask]), color=mineos_plot_colour, linewidth=lw_syn, label='MINEOS', alpha=alpha3Dcodes)
    if ictr == 0:
        leg_lines.append(lnm)
        leg_lines_3D.append(lnm3d)

    # Plot each code for this station
    for icode in [axitra, qssp, spec]:

        ivarax = 0
        for ivar in [icode.z_acc[stn], icode.z_grav[stn], icode.z_pegs[stn]]:
            ltmp, = ax[ivarax].plot(icode.time[stn], 9 + ivar - ictr, color=icode.plot_colour, linewidth=lw_syn, label=icode.plot_label)
            ivarax += 1

        # Only add to legend string once (not for every station)
        if ictr == 0:
            leg_lines.append(ltmp)


    # For comparison with 3D S40RTS:
    iic = 0
    for icode in [axitra, qssp, spec, spec3D]:
        # S40RTS needs time cut
        tmask = icode.time[stn] > 0

        # Alter the alpha value depending on the code
        if icode.plot_label != 'S40RTS':
            ialpha = alpha3Dcodes
            lsty ="-"
        else:
            ialpha = 1
            lsty='-'

        ltmp3D, = ax3D.plot(icode.time[stn][tmask], 9 + icode.z_pegs[stn][tmask] - ictr, color=icode.plot_colour,
                            linewidth=lw_syn, linestyle=lsty, label=icode.plot_label, alpha=ialpha)
        # Only add to legend string once (not for every station)
        if ictr == 0:
            leg_lines_3D.append(ltmp3D)
        iic +=1

    # Add station name and horizontal lines
    for a in [ax3D, ax[0]]:
        a.text(x=-30, y= 9 - ictr, s=stn, weight="bold", horizontalalignment='center', fontsize=10)
        a.text(x=-30, y= 9 - ictr-0.2, s=f"{int(epidists[stn])} km", weight="bold", horizontalalignment='center', fontsize=7.5)
        #a.axhline(y=9 - ictr, xmin=0/450, xmax=100/450, color='grey', linewidth=1, alpha=0.3)



    ictr += 1



# Format legend
leg_titles    = ['MINEOS', 'AXITRA', 'QSSP', 'SPECFEM', 'Real data']
leg_titles_3D = ['MINEOS', 'AXITRA', 'QSSP', 'SPECFEM', 'S40RTS', 'Real data']

leg_lines.append(lreal)
leg_lines_3D.append(lreal)

legend   = ax[-1].legend(leg_lines,leg_titles, fontsize=8)
legend3D = ax3D.legend(leg_lines_3D,leg_titles_3D,  fontsize=8)
for l in [legend, legend3D]:
    frame = l.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')
    l.get_frame().set_alpha(None)

# Formatting the limits, labels etc
for a in all_axes:
    a.set_ylim([0, 9.5])
    a.set_xlim([0, 350])
    a.axvline(0, ymin=vlinmin, ymax=1, linewidth=2, color='k')
    a.set_xlabel('Time since rupture initiation [s]', weight='bold')

    # Remove standard labels
    for xlabel_i in a.axes.get_yticklabels():
        xlabel_i.set_fontsize(10.0)
        xlabel_i.set_visible(False)

    # Scale bar for plot
    a.errorbar(x=300, y=2.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
    a.text(x=305, y=2.4, s=r'1 $nm/s^2$' )

    # Remove ugly spines!
    a.tick_params(labelleft=False, left=False)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.spines['left'].set_visible(False)

    # X tick labels
    a.set_xticks([0, 100, 200, 300])
    a.set_xticklabels(('0', '100', '200', '300'))



ax[0].set_title('Gravitational acceleration', weight='bold')
ax[1].set_title('Elastic acceleration', weight='bold')
ax[2].set_title('PEGS signal', weight='bold')



subplotabc = 'abc'
for i in range(3):
    ax[i].text(x=5 , y=vlinmin*9.5, s=f"({subplotabc[i]})", weight='bold', fontsize=10)


# Save pdf and output to screen
for suffix in ['png', 'pdf']:
    path = f'../../Figures/{suffix}s/'
    fig.savefig(f'{path}/PEGS_tohoku_benchmark.png', format=suffix)
    fig3D.savefig(f'{path}/PEGS_tohoku_3D_compare.png', format=suffix)
plt.show()

