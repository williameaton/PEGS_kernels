import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec
from wetools.plotting import setup_we_mpl
setup_we_mpl()


prem1c = np.loadtxt('prem_1_crust.txt')
prem2c = np.loadtxt('prem_2_crust.txt')
ak135w = np.loadtxt('ak135_wei.txt')


figgrid = plt.figure(constrained_layout=True, figsize=(15,7))
spec = gridspec.GridSpec(ncols=10, nrows=15, figure=figgrid)
axfull = figgrid.add_subplot(spec[:, :2])

axes = []
for row  in range(3):
    rowax = []
    for cols in range(4):
        rowax.append(figgrid.add_subplot(spec[(row*5)  : (row+1)*5 , 2+ 2*cols : 4+ 2*cols  ]))
    axes.append(rowax)
axes = np.array(axes)
#fig, ax = plt.subplots(1, 3, sharey=True)
#figrho, axrho = plt.subplots(1, 4, sharey=True)


pltclrs = ['#004488', '#DDAA33', '#BB5566']
mctr = 0
for m in [prem1c, prem2c, ak135w]:
    for iax in range(3):
        axfull.plot(m[:,2+iax], m[:,1], color=pltclrs[mctr])
        axfull.set_ylim([1000, 0])

    mctr += 1
axfull.set_ylabel('Depth [km]')
axfull.set_xlabel(r'Velocity [$km\,s^{-1}$] / Density [$g\, cm^{-3}$]')

# This was originally written for rho but now loops for vs and vp too:


for ivar in range(3):
    # I want to piecewise interpolate each of these
    interp_models = []
    for m in [prem1c, prem2c, ak135w]:

        mask = m[:, 1] <= 400 # top 250 km

        depth = m[mask, 1]
        rho   = m[mask, ivar+2]

        breaks = np.where(np.diff(depth) <= 0)[0] + 1
        segments = np.split(np.arange(len(depth)), breaks)

        npts = 10009

        new_depth = np.linspace(0, 250, npts)
        new_rho   = np.full_like(new_depth, np.nan)

        # Loop through each segment and interpolate
        for seg in segments:
            dseg = depth[seg]
            vseg = rho[seg]

            # Avoid issues with repeated depths in segment
            unique_d, idx = np.unique(dseg, return_index=True)
            unique_v = vseg[idx]

            # Create interpolator
            interp = interp1d(unique_d, unique_v, kind='linear', bounds_error=False)

            # Find points in new_depth that fall within the current segment
            mask = (new_depth >= unique_d[0]) & (new_depth <= unique_d[-1])
            new_rho[mask] = interp(new_depth[mask])

        interp_models.append(new_rho)


    # axes for this one:
    axrho = axes[ivar, :]

    mctr = 0
    for i in range(3):
        axrho[0].plot(interp_models[i], new_depth, '-', color=pltclrs[i])
        axrho[1+i].axvline(0, alpha=0.5, color='k', linestyle ='--')
        axrho[1+i].axhline(20, alpha=0.5, color='k', linestyle =':')

        axrho[1 + i].set_xlim([-25, 25])
        axrho[1 + i].set_ylim([250, 0])
        axrho[1 + i].set_yticklabels([])

        if ivar==2:
            axrho[1 + i].set_xlabel('% model difference')

        # Add a shaded patch for the 15- 24.4 km region:
        poly_coords = [
            (-25, 15.0), (25, 15.0),
             (25, 24.4), (-25, 24.4)
        ]
        axrho[1+i].add_patch(
            plt.Polygon(poly_coords, color='#DDAA33',
                        alpha=0.7)
        )

    axrho[0].set_ylim([250, 0])
    axrho[0].legend(['PREM 1C', 'PREM 2C', 'AK135'])

    axrho[1].plot(100*(interp_models[1]-interp_models[0])/interp_models[0], new_depth, color='k', linestyle='-')
    #axrho[1].plot((interp_models[1]-interp_models[0])/interp_models[0], new_depth, color=pltclrs[1], linestyle='--')

    axrho[2].plot(100*(interp_models[2]-interp_models[0])/interp_models[0], new_depth, color='k', linestyle='-')
    #axrho[2].plot((interp_models[2]-interp_models[0])/interp_models[0], new_depth, color=pltclrs[2], linestyle='--')


    axrho[3].plot(100*(interp_models[2]-interp_models[1])/interp_models[1], new_depth, color='k', linestyle='-')
    #axrho[3].plot((interp_models[2]-interp_models[1])/interp_models[1], new_depth, color=pltclrs[2], linestyle='--')

labels = ['PREM 1C\nvs\nPREM 2C', 'PREM 1C\nvs\nAK135+WEI', 'PREM 2C\nvs\nAK135+WEI']

for j in range(3):
    for i in range(3):
        axes[i,j+1].text(x =-15, y=235, s=labels[j], fontsize=8, horizontalalignment='center')

for row in range(3):
    axes[row,0].set_ylabel('Depth [km]')

# Add row text:
fs = 14
xx = 14
axfull.text(x=xx, y=175, s=r'$V_p$', fontsize=fs)
axfull.text(x=xx, y=520, s=r'$V_s$', fontsize=fs)
axfull.text(x=xx, y=870, s=r'$\rho$', fontsize=fs)


figgrid.savefig('models.pdf', format='pdf')

plt.show()