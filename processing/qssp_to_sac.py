import obspy
import numpy as np
from obspy.core.trace import Trace

qssp_dir = '/Users/eaton/Documents/Princeton/PEGS_kernels/forward_simulations/Tohoku/data/qssp/prem_elastic/'
raw_dir  = qssp_dir + '/raw'
grav_dir = qssp_dir + '/Grav'
zacc_dir = qssp_dir + '/Zacc'
name     = 'ps_vapf_0d'

trr = Trace()

chl = 'z'

# Loop for grav and ground acc
for dataloop in range(2):

    # Data file name:
    if dataloop==0:
        tmp_fname = f'{raw_dir}/{name}_grav_acce_{chl}.dat'
        outdir = grav_dir
    else:
        tmp_fname = f'{raw_dir}/{name}_acce_{chl}.dat'
        outdir = zacc_dir

    # Read station order
    stns = open(tmp_fname, 'r').readline().split()[1:]

    # Load data:
    data = np.loadtxt(tmp_fname, skiprows=1)

    # Check loaded data for each station
    nstns = len(stns)
    assert(np.shape(data)[1] == nstns + 1)

    # Set time data
    time = data[:,0]
    dt = np.mean(time[1:] - time[:-1])


    # Loop each station
    for istn in range(nstns):
        trr.data = data[:, istn+1]
        trr.stats.delta = dt

        trr.write(f"{outdir}/{stns[istn]}.sac", format="SAC")


