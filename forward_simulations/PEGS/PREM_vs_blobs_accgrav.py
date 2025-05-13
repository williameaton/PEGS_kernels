# Plots Tohoku PEGS from SPECFEM and QSSP compared to tests from the 'blobs' directory
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../classes')
from pegs import create_pegs
import pickle

# Some plotting stuff:
leg_titles = ['PREM 2C Mw 8.5',
              'SPECFEM PREM 2C',
              'QSSP PREM 2C',
              'SPECFEM PREM 1C',
              'QSSP PREM 1C',
              ]

# Data directories
ddir      = '../../data/data_forward_pegs'
blobs_dir = '../../data/data_blobs/'
qssp_dir  = f'{ddir}/qssp'

figsize  = (13, 7)         # Figure size
lw_syn   = 2               # Synthetic data line width
vlinmin  = 0.015           # Min y value of vertical lines

# Load the epicentral distances in km for annotation:
with open('../../Figures/data_for_plots/epi_dist.pkl', 'rb') as fp:
    epidists = pickle.load(fp)

# Colourscheme
medcontrast = ['#EECC66', '#EE99AA', '#6699CC',
               '#997700', '#994455', '#004488']



# PREM 1 crust
prem_1c = create_pegs(f'{ddir}/PREM_1C', code='SPECFEM')
prem_1c.load()
prem_1c.plot_colour = medcontrast[5]

# PREM 2 crust
prem_2c = create_pegs(f'{ddir}/PREM_2C', code='SPECFEM')
prem_2c.load()
prem_2c.plot_colour = medcontrast[4]



# Load elastic prem with physical dispersion
qssp_PREM_1c = create_pegs(f'{qssp_dir}/prem_elastic_1C/', code='QSSP')
qssp_PREM_1c.load()
qssp_PREM_1c.cut_to_plot()
qssp_PREM_1c.plot_colour = medcontrast[2]


qssp_PREM_2c = create_pegs(f'{qssp_dir}/prem_elastic_2C/', code='QSSP')
qssp_PREM_2c.load()
qssp_PREM_2c.cut_to_plot()
qssp_PREM_2c.plot_colour = medcontrast[1]


# Wrong magnitude - 8.5
qssp_wrongmag = create_pegs(f'{qssp_dir}/prem_elastic_2C_Mag8.5/', code='QSSP')
qssp_wrongmag.load()
qssp_wrongmag.cut_to_plot()
qssp_wrongmag.plot_colour = medcontrast[0]



# AK135 anelastic no physical dispersion
# qssp_ak135_el = create_pegs(f'{qssp_dir}/original_AK135/elastic_output/', code='QSSP')
# qssp_ak135_el.load()
# qssp_ak135_el.cut_to_plot()
# qssp_ak135_el.plot_colour = medcontrast[1]

# # SPECFEM AK135-WEI
# ak135wei = create_pegs(f'{ddir}/ak135wei', code='SPECFEM')
# ak135wei.load()
# ak135wei.plot_colour = medcontrast[4]


# # Test 5 - altering vp/vs/rho for central 450 km cylinder
# cylinder_mp_vpvs = create_pegs(f'{blobs_dir}/test_5_central_cylinder/forward_altervpvsrho/', code='SPECFEM')
# cylinder_mp_vpvs.load()
# cylinder_mp_vpvs.plot_colour = '#DDAA33'
#
# # Test 6 - 75 km at source
# cylindersrc_vpvs = create_pegs(f'{blobs_dir}/test_6_source_cylinder/', code='SPECFEM')
# cylindersrc_vpvs.load()
# cylindersrc_vpvs.plot_colour = '#004488'


simulations = [qssp_wrongmag,
               prem_2c, qssp_PREM_2c,
               prem_1c, qssp_PREM_1c
               ]

# Create figures and set tight layouts:
fig, axes     = plt.subplots(1,3, figsize=figsize)
fig.set_tight_layout(True)

# Stores line objects for mpl legend in loop
leg_lines    = []

# Loop over the grav and zacc:
for iax in range(3):
    if iax == 0:
        plot_grav = True
        plot_pegs = False
        axtitle = 'Gravitational acceleration'
    elif iax == 1:
        plot_grav = False
        plot_pegs = False
        axtitle = 'Elastic acceleration'
    else:
        plot_grav = False
        plot_pegs = True
        axtitle = 'PEGS signal'
    ax = axes[iax]

    # Loop over each station
    ictr = 0
    for stn in simulations[0].stns:

        # Plot each code for this station
        for icode in simulations:

            if plot_grav:
                plotvar = icode.z_grav[stn]
            else:
                plotvar = icode.z_acc[stn]

            if plot_pegs:
                plotvar = icode.z_pegs[stn]


            ltmp, = ax.plot(icode.time[stn], 9 + plotvar - ictr, color=icode.plot_colour, linewidth=lw_syn)
            # Only add to legend string once (not for every station)
            if ictr == 0:
                leg_lines.append(ltmp)

        # For left-most axis:
        if iax == 0:
            # Add station name and horizontal lines
            ax.text(x=-50, y= 9 - ictr, s=stn, weight="bold", horizontalalignment='center', fontsize=10)
            ax.text(x=-50, y= 9 - ictr-0.2, s=f"{int(epidists[stn])} km", weight="bold", horizontalalignment='center', fontsize=7.5)


        ictr += 1
    ax.set_title(axtitle, weight='bold')

    # Format legend
    legend   = ax.legend(leg_lines, leg_titles, fontsize=8)
    for l in [legend]:
        frame = l.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')
        l.get_frame().set_alpha(None)

    # Formatting the limits, labels etc
    for a in [ax]:
        a.set_ylim([0, 9.5])
        a.set_xlim([0, 350])
        a.axvline(0, ymin=vlinmin, ymax=1, linewidth=2, color='k')
        a.set_xlabel('Time since rupture initiation [s]', weight='bold')

        # Remove standard labels
        for xlabel_i in a.axes.get_yticklabels():
            xlabel_i.set_fontsize(10.0)
            xlabel_i.set_visible(False)

        # Scale bar for plot
        a.errorbar(x=300, y=3.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
        a.text(x=305, y=3.4, s=r'1 $nm/s^2$' )

        # Remove ugly spines!
        a.tick_params(labelleft=False, left=False)
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.spines['left'].set_visible(False)

        # X tick labels
        a.set_xticks([0, 100, 200, 300])
        a.set_xticklabels(('0', '100', '200', '300'))

plt.savefig('../../Figures/pdfs/elastic_grav_acc_codes.pdf')

plt.show()