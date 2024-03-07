# Plots a histogram of the magnitude of values for a kernel, globally

import ensightreader
import matplotlib.pyplot as plt
import numpy as np

fpath ="/Users/eaton/Downloads/download/download/"
varname = "grav1_kernel"


dpos_tot  = np.array([])
dneg_tot  = np.array([])
dfull_tot = np.array([])


for proc in range(864):

    case = ensightreader.read_case(f"{fpath}/reg1_proc{proc}.case")
    geofile = case.get_geometry_model()

    part_names = geofile.get_part_names()
    part = geofile.get_part_by_name(part_names[0])
    N = part.number_of_nodes

    #with geofile.open() as fp_geo:
    #    node_coordinates = part.read_coordinates(fp_geo)  # np.ndarray((N, 3), dtype=np.float32)

    variable = case.get_variable(varname)


    with open(variable.file_path, "rb") as fp_var:
        d = variable.read_node_data(fp=fp_var, part_id=part.part_id)

    # split into positive and negative:
    dpos = np.log10(d[d>=0])
    dneg = np.log10(np.abs(d[d<0]))

    dtot = np.concatenate((dpos, dneg))

    dpos_tot = np.concatenate((dpos_tot, dpos))
    dneg_tot = np.concatenate((dneg_tot, dneg))
    dfull_tot = np.concatenate((dfull_tot, dtot))



maxval = np.amax(dfull_tot)


fig, ax = plt.subplots(1,2)

upper = np.ceil(maxval)
bins = np.arange(-23, upper, 0.5)
ax[0].hist(dpos_tot, bins=bins, alpha=0.5, histtype='step', color='red')
ax[0].hist(dneg_tot, bins=bins, alpha=0.5, histtype='step', color='blue')
ax[0].hist(dfull_tot, bins=bins, alpha=1, histtype='step', color='black')

ax[0].set_title(f'Processor {proc}: Variable {varname}')
ax[0].set_xlabel('Log amplitude of corner GLL point')
ax[0].set_ylabel('Count')
ax[0].set_xlim([-23, upper -0.5])

ax[0].axvline(maxval, color='red') # max val
ax[1].axvline(maxval, color='red') # max val


# remove anything that is less than some threshold orders of mag below the max value:
threshold = 3
thr = maxval - threshold

dfull_thresh = dfull_tot[dfull_tot > thr]

ax[0].axvline(thr, color='teal')


ax[1].hist(dfull_thresh, bins=int(1+np.log2(len(dfull_thresh))), alpha=1, histtype='step', color='green')


ax[1].set_xlim([np.min(dfull_thresh)-0.5, np.max(dfull_thresh)+0.5])


plt.show()