# loads kernel values from MDJ ensight files for pegs kernel and plots the value decay with distance
import ensightreader
import matplotlib.pyplot as plt
import numpy as np
from wetools.plotting import setup_we_mpl
setup_we_mpl()

src_loc = [-0.6338467052282285, 0.47677109139887874, 0.6090383244736256]

procs = np.loadtxt('ensight/wanted_procs.txt').astype(int)

fig, ax = plt.subplots(2,3, figsize=(13,7), sharex=True, sharey=True)

knls = ["kappa", "mu", "rhonotprime", "alpha", "beta", "rho"]
knltitle = [r'$\kappa$', r'$\mu$', r'$\rho$', r'$\alpha$', r'$\beta$', r"$\rho'$"]

icol = 0
irow = 0
for iknl in range(6):

    knlname = knls[iknl]

    radius = []
    kernel = []
    for iproc in procs:
        casefile = f'ensight/reg1_proc{iproc}.case'
        case = ensightreader.read_case(casefile)
        geofile = case.get_geometry_model()

        part_names = geofile.get_part_names()
        part = geofile.get_part_by_name(part_names[0])
        N = part.number_of_nodes

        with geofile.open() as fp_geo:
            node_coordinates = part.read_nodes(fp_geo)

        # Read data for this kernel
        variable = case.get_variable(knlname)
        with open(variable.file_path, "rb") as fp_var:
            vardata = variable.read_node_data(fp_var, 1)

        # We now need to compute the distance from the source in cartesian coordinates

        dist = ((node_coordinates[:,0] - src_loc[0])**2 +
                (node_coordinates[:,1] - src_loc[1])**2 +
                (node_coordinates[:,2] - src_loc[2])**2)**0.5

        mask = dist < 250/6371
        if len(dist[mask]) > 0:
            radius = radius +  list(dist[mask][:])
            kernel = kernel + list(vardata[mask][:])

    radius = np.array(radius)
    kernel = np.array(kernel)

    knlplot=np.log10(np.abs(kernel))

    ax[irow, icol].scatter(radius, knlplot, alpha=0.05)

    ax[irow, icol].set_title(knltitle[iknl], weight='bold')


    ind = 3
    rr = np.linspace(0.00001, 0.04, 1000)
    yy = np.log10((1/(rr**ind)))
    yy = yy + np.max(knlplot) - 9
    ax[irow, icol].plot(rr,yy, 'r')

    ind = 4
    rr = np.linspace(0.00001, 0.04, 1000)
    yy = np.log10((1/(rr**ind)))
    yy = yy + np.max(knlplot) - 11
    ax[irow, icol].plot(rr,yy, 'orange')


    ind = 5
    rr = np.linspace(0.00001, 0.04, 1000)
    yy = np.log10((1/(rr**ind)))
    yy = yy + np.max(knlplot) - 13
    ax[irow, icol].plot(rr,yy, 'purple')

    icol += 1

    if icol==3:
        irow += 1
        icol = 0


ax[0, 0].set_ylim([-28, -18])

dist_ticks = np.arange(0, 275, 50)
for icol in range(3):
    ax[-1, icol].set_xticks(dist_ticks/6371)
    ax[-1, icol].set_xticklabels(dist_ticks)
    ax[-1, icol].set_xlabel('Distance from source [km]')

for irow in range(2):
    ax[irow, 0].set_ylabel('Log Absolute kernel magnitude')


fig.savefig(f"knl_decay.pdf", format='pdf')
plt.show()