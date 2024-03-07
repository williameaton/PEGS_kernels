import matplotlib.pyplot as plt
import numpy as np
from wetools import norm


glb_dir = './globe_forward_output'
glb_adj_dir = './globe_adjoint_output'
spfmx_dir = './spfmx_forward_output'
spfmx_adj_dir = './spfmx_adjoint_output'
chls = ['Z', 'N', 'E']
depth = '0'
stn = 'A4'


clr_sp     = 'slategrey'
clr_sp_adj = 'k'
clr_glb    = 'darkorchid'
clr_glb_adj= 'coral'


fig = plt.figure(constrained_layout=True, figsize=(12, 7))
gs = fig.add_gridspec(3, 3)


# Source time function axis:
ax_stf = fig.add_subplot(gs[:2, 2])
tshift = 50
ax_stf.axvline(0, color='k', linestyle='--', alpha=0.5)
# Load
glb_stf = np.loadtxt(f"{glb_dir}/plot_source_time_function.txt")
spfmx_stf = np.loadtxt(f"{spfmx_dir}/plot_source_time_function.txt")

from copy import copy
backstf = copy(glb_stf)
backstf[:, 1] = backstf[::-1, 1]

# Plot
ax_stf.plot(glb_stf[:,0], glb_stf[:,1]*glb_stf[:,2], color=clr_glb)
ax_stf.plot(backstf[:,0], backstf[:,1]*backstf[:,2], color=clr_glb_adj)
#ax_stf.plot(spfmx_stf[:,0] - tshift, spfmx_stf[:,1], 'r--')
ax_stf.set_xlabel("Time")


# Plot traces:
axZ = fig.add_subplot(gs[0, :2])
axN = fig.add_subplot(gs[1, :2])
axE = fig.add_subplot(gs[2, :2])
axchls = [axZ, axN, axE]

nochls = len(chls)
for ichl in range(nochls):
    chl = chls[ichl]

    # Load:
    glob  = np.loadtxt(f"{glb_dir}/IU.{stn}_{depth}.MX{chl}.sem.ascii")
    spfmx = np.loadtxt(f"{spfmx_dir}/{stn}_{depth}.IU.MX{chl}D.sem.ascii")

    lg, = axchls[ichl].plot(glob[:,0],  norm(glob[:,1]), color=clr_glb)
    ls, = axchls[ichl].plot(spfmx[:,0], norm(spfmx[:,1])-0.1, color=clr_sp)

    glob_adj = np.loadtxt(f"{glb_adj_dir}/IU.{stn}_{depth}.MX{chl}.sem.ascii")
    spfmx_adj = np.loadtxt(f"{spfmx_adj_dir}/{stn}_{depth}.IU.MX{chl}D.sem.ascii")

    lga, = axchls[ichl].plot(glob_adj[:,0], norm(glob_adj[:,1]) + 0.1, color=clr_glb_adj)
    lgs, = axchls[ichl].plot(spfmx_adj[:,0], norm(spfmx_adj[:,1]) + 0.2, color=clr_sp_adj)

    axchls[ichl].set_title(f'Channel: {chl}')

    axchls[ichl].legend([lg, ls, lga, lgs],['Globe forward', 'SPFMX forward', 'Globe tricked adjoint', 'SPFMX tricked adjoint'], loc='upper left')

# Plot gravitational potential trace for both
ax_G = fig.add_subplot(gs[2, 2])
spfmx = np.loadtxt(f"{spfmx_dir}/{stn}_{depth}.IU.MXG.sem.ascii")
spfmx_adj = np.loadtxt(f"{spfmx_adj_dir}/{stn}_{depth}.IU.MXG.sem.ascii")
ax_G.plot(spfmx[:,0], norm(spfmx[:,1]), color=clr_sp)
ax_G.plot(spfmx_adj[:,0], norm(spfmx_adj[:,1]), color=clr_sp_adj, linestyle=':')
ax_G.legend([r'Forward $\phi$', r'Tricked adjoint $\phi$'])
plt.show()