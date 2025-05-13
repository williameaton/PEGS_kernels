import sys
from copy import copy
import numpy as np
sys.path.append(f"../../../classes")
sys.path.append(f"../../../Figures")
from finite_difference import FD, FdStation
from plot_all_stations import plot_all_stations
from wetools.plotting import save_figs_to_single_pdf, setup_we_mpl, subplots_hide_xaxes, align_y_labels
import matplotlib.pyplot as plt
from colour_schemes import hex

share_stn_ylims = False

alpha_for_alt_gradients = 1
xticklabelsize = 7


# Setup my standard MPL
setup_we_mpl()

# Event location and channels recorded
src_lat  = 37.5200
src_lon  = 143.05
acc_chls = 'ZNE'

zhang_ddir = '../../data/ZHANG_data/Data_Zhang2023/WEHS/'


# Multiplier value for acceleration
mult = -1
scalepower = 15
scaling = 10**scalepower

# Plotting colours
HC = hex.Hexes().hex['HighContrast_PT']

nplotchls = 5

figs = []; figzz = []

for NEX in [384]:
    for dx in [5000]:

        spfm_data_path = f'./data/NEX{NEX}_w_perfsphere/outputs/'

        # Create a Finite Diff setup from the README parameters
        fd = FD(from_readme=True, readme_fname='qssp_zhang_stns_plus_PEGSstns/README')
        fd.reorder_stations_by_distance(src_lat, src_lon)

        #fd.save_distances_to_pickle('../../../Figures/data_for_plots/epi_dist.pkl')

        # Figures
        fig, ax = plt.subplots(fd.nstns, nplotchls, figsize=(16,9.5))
        fig.suptitle(f"Gradient of gravitational acceleration (NEX {NEX})", fontsize=15, weight='bold', y=0.98)


        # Load grav acceleration and potential data
        fd.load_specfem_data(types=['pgrav', 'gpot'], dpath=spfm_data_path,
                             acc_mult=mult)
        fd.load_qssp_grav_acc(fpath=zhang_ddir)


        # -------------------------------- QSSP STRAIN --------------------------------
        # Compute strain from centred FD3 of QSSP acc
        chls =  ['NN', 'EE', 'ZN', 'EN', 'ZE']
        fd.strain_from_acc(method=f'centred_FD3', dx=dx, code='qssp')
        axi = 0
        for istn in fd.stations:
            ichl = 0
            d = istn.qssp
            time = d['time']
            for chl in chls:
                lqssp, = ax[axi, ichl].plot(time, scaling *  d[f'{chl}.STRAIN'])
                ichl += 1
            axi += 1

        # -------------------------------- SPECFEM STRAIN --------------------------------
        # Compute strain from centred FD5 of acc
        chls =  ['NN', 'EE', 'ZN', 'EN', 'ZE']
        fd.strain_from_acc(method=f'centred_FD5', dx=dx)
        axi = 0
        for istn in fd.stations:
            ichl = 0
            d = istn.data
            time = d['time']
            for chl in chls:
                lfd5, = ax[axi, ichl].plot(time, scaling*d[f'{chl}.STRAIN'], alpha=alpha_for_alt_gradients)
                ichl +=1
            axi += 1

        # Compute strain from 2nd order FD of gpot:
        # Note that only NN and EE are being updated here and the other channels are just from the computation above
        fd.strain_from_gpot(method=f'centred_FD5', dx=dx)
        axi = 0
        for istn in fd.stations:
            ichl = 0
            d = istn.data
            time = d['time']
            for chl in chls:

                if ichl > 1:
                    mmult = 1
                else:
                    mmult = mult

                # We need it to plot to be consistent with the legend colours but want it to be a null plot (not visible
                # Hence the alpha = 0
                # the values for channels other than NN and EE are inherited from previous method above
                if chl == 'NN' or chl == 'EE':
                    lpot, = ax[axi, ichl].plot(time,  mmult * scaling * d[f'{chl}.STRAIN'], alpha=alpha_for_alt_gradients)
                else:
                    ax[axi, ichl].plot(time,  0* mmult * scaling * d[f'{chl}.STRAIN'], alpha=0)


                ichl +=1
            axi += 1

        # Load strain directly:
        chls =  ['NN', 'EE', 'ZN', 'EN', 'ZE']
        fd.load_specfem_data(types=['strain'], dpath=spfm_data_path)
        ichl = 0
        for chl in chls:
            axi = 0
            for istn in fd.stations:
                d = istn.data
                time = d['time']

                ldirect, = ax[axi, ichl].plot(time, mult * scaling * d[istn.name + '_' + str(istn.cstn)][f'{chl[::-1]}.STRAIN'])
                axi += 1
            ichl +=1
        # -------------------------------- END OF SPECFEM STRAIN --------------------------------

        for i in range(fd.nstns):
            ax[i, 0].set_ylabel(f'{fd.stations[i].name}\n{int(fd.stations[i].epi_dist_km)} km', weight='bold', fontsize=9, rotation=0)
            ax[i, 0].yaxis.set_label_coords(-0.35, 0.2)

        subplots_hide_xaxes(ax, nrows=fd.nstns, ncols=nplotchls, keep_lowest=False)


        chls =  ['NN', 'EE', 'ZN', 'EN', 'ZE']
        for j in range(nplotchls):
            ax[0, j].set_title(f'Channel {chls[j]}', weight='bold')


        # Set legend:
        leglines = [lqssp, lfd5, lpot, ldirect]
        leglbls  = ['QSSP', '1st order FD5 of acceleration', '2nd order FD5 of potential', 'SPECFEM output']
        nlines   = len(leglines)
        assert(len(leglbls)==nlines)

        ax[-1,2].legend(leglines, leglbls,
                        ncol=nlines,
                        loc='lower center',
                        bbox_to_anchor=(0.5, -1.25))

        for istn in range(fd.nstns):
            for ichl in range(nplotchls):
                ax[istn, ichl].set_xlim([0, fd.stations[istn].time_cutoff])

                ax[istn, ichl].tick_params(axis='y', labelsize=xticklabelsize)

            if share_stn_ylims:
                minlim = 100
                maxlim = -100
                for axi in ax[istn,:]:
                    ylim = axi.get_ylim()

                    ymin = ylim[0]
                    ymax = ylim[1]

                    if ymin < minlim:
                        minlim = ymin
                    if ymax > maxlim:
                        maxlim = ymax
                # Set the final value:
                for axi in ax[istn, :]:
                    axi.set_ylim([minlim, maxlim])


        ax[-1, 0].set_xlabel(f'All y-axes scaled by $10^{{{int(scalepower)}}}$', fontsize=8)
        ax[-1, 0].xaxis.set_label_coords(0.35, -0.25)

        figs.append(fig)



plt.show()

#
# map = plot_all_stations('../../../Figures/data_for_plots/PEGS_STATIONS')
#
# finalfigs = [map]
# for j in [figs, figzz]:
#     for i in j:
#         finalfigs.append(i)
# save_figs_to_single_pdf(finalfigs, pdf_name='compare_FDs_reprod.pdf')
#
#
# # Since the NEX 384 is the last figure we can save it separately as this
# fig.subplots_adjust(left=0.1, right=0.99, top=0.91, bottom=0.05, hspace=0.4, wspace=0.4)
#
fig.savefig('NEX384_strain_benchmark.pdf', format='pdf')


