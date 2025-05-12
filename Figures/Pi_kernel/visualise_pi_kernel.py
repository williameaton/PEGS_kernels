import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from wetools.plotting import setup_we_mpl
setup_we_mpl()
def norm(x):
    return x/np.max(np.abs(x))

# Specify a receiver based on its cartesian coordinates:
# MD
xrcvr = -0.453641593584130
yrcvr  = 0.548537798585614
zrcvr  = 0.702364284465902
stnnorm = np.array([xrcvr, yrcvr, zrcvr])

# Create a sphere
r = 1
pi = np.pi
cos = np.cos
sin = np.sin
n = 100
phi, theta = np.mgrid[0:np.pi/2:n*1j, 0:2*np.pi:n*1j]
x = r*sin(phi)*cos(theta)
y = r*sin(phi)*sin(theta)
z = r*cos(phi)




# For each point on the x,y,z compute the Pi matrix:
PI = np.zeros((3,3,n, n))
TotalPiMatrix = np.zeros((3,3,n, n))
xr = x - xrcvr
yr = y - yrcvr
zr = z - zrcvr
dr = [xr, yr, zr]

rr = (xr**2 + yr**2 + zr**2)**0.5
ir3 = 1/rr**3
ir5 = 1/rr**5


for i in range(3):
    for j in range(3):
        PI[i,j,:,:]            =  3*dr[i] * dr[j] * ir5
        TotalPiMatrix[i,j,:,:] = -3*dr[i] * dr[j] * ir5

    TotalPiMatrix[i,i,:,:] += ir3

# For colormap of the log, we want to round to the nearest whole? number
"""for mmm in range(3):
    for nnn in range(3):
        if (mmm == nnn):
            crosslog10 = np.log10(PI[mmm,nnn,:,:])
            Iplot      = np.log10(ir3)
            TotalPi    = ir3 - PI[mmm,nnn,:,:]
        else:
            Iplot = 0*ir3
            crosslog10 = 0*ir3
            TotalPi = - PI[mmm,nnn,:,:]

        # To compute total Pi we need to mask out the negative values so we can take a log
        LogTotalPi = TotalPi*0
        LogTotalPi[:,:] = np.NaN
        TotalPi_mask = TotalPi > 0
        LogTotalPi[TotalPi_mask] = np.log10(TotalPi[TotalPi_mask])


        # Similarly we now want to get the negative values but make them positive for plotting:
        NegLogTotalPi = TotalPi*0 + np.NaN
        NegPi_mask    = TotalPi < 0
        NegLogTotalPi[NegPi_mask] = np.log10(-TotalPi[NegPi_mask])

        maxLogTotal_nonan    = np.nanmax(LogTotalPi)
        maxNegLogTotal_nonan = np.nanmax(NegLogTotalPi)


        # Get the min and max of the two plots:

        vmax3 = np.floor(np.max([np.max(Iplot), np.max(crosslog10), maxLogTotal_nonan, maxNegLogTotal_nonan]))
        vmin3 = vmax3 - 8 #np.floor(np.min(log10))
        vspan = vmax3 - vmin3

        print('Vmax: ', vmax3)

        # Scale each plot
        crosslog10 = (crosslog10 - vmin3)/vspan
        Iplot = (Iplot - vmin3)/vspan
        LogTotalPi = (LogTotalPi - vmin3)/vspan
        NegLogTotalPi = (NegLogTotalPi - vmin3)/vspan



        vrange = np.arange(vmin3, vmax3+1)
        nints = len(vrange)

        # Create custom colorbar
        ni = nints-1  # Number of discrete colors you want
        oranges = plt.cm.Oranges(np.linspace(0, 1, ni))  # Sample n colors from a continuous colormap
        blues   = plt.cm.Blues(np.linspace(0, 1, ni))  # Sample n colors from a continuous colormap
        discrete_orange = LinearSegmentedColormap.from_list("orange_disc", oranges, N=ni)
        discrete_blue   = LinearSegmentedColormap.from_list("orange_disc", blues, N=ni)


        fig = plt.figure()
        G = gridspec.GridSpec(8,6)

        cbax_or = fig.add_subplot(G[-1,4:])
        cbax_bl = fig.add_subplot(G[-1,:2])

        ax_I        = fig.add_subplot(G[1:7,:2], projection='3d')
        ax_cross    = fig.add_subplot(G[1:7,2:4], projection='3d')
        ax_complete = fig.add_subplot(G[1:7,4:6], projection='3d')




        IplotMap    = ax_I.plot_surface(x, y, z,        facecolors=discrete_orange(Iplot))
        CrossMap    = ax_cross.plot_surface(x, y, z,    facecolors=discrete_blue(crosslog10))

        CompleteMap = ax_complete.plot_surface(x, y, z, facecolors=discrete_orange(LogTotalPi))
        CompleteMap = ax_complete.plot_surface(x, y, z, facecolors=discrete_blue(NegLogTotalPi))


        cbar_range = np.arange(vmin3, vmax3+1)
        mappable = plt.cm.ScalarMappable(cmap=discrete_orange)
        mappable.set_array(cbar_range)
        mappable.set_clim(0, 1)  # Normalized range [0, 1]
        fig.colorbar(mappable, cax=cbax_or, orientation='horizontal')

        mappable = plt.cm.ScalarMappable(cmap=discrete_blue)
        mappable.set_array(cbar_range)
        mappable.set_clim(0, 1)  # Normalized range [0, 1]
        fig.colorbar(mappable, cax=cbax_bl, orientation='horizontal')
        cbax_bl.set_xlim([1,0])

        cbar_ticks = (cbar_range - vmin3)/vspan
        cbax_or.set_xticks(cbar_ticks, cbar_range.astype(int))
        cbax_bl.set_xticks(cbar_ticks, cbar_range.astype(int)[::1])

        cbax_or.set_xlabel('Positive-value log10 amplitude')
        cbax_or.set_xlabel('Negative-value log10 amplitude')

        ax_I.set_title(r"$||\mathbf{r} - \mathbf{r}'||^{-3}$")

        ax_cross.set_title(
            rf"$ \frac{{ -3 (r_{{{mmm}}} - r_{{{mmm}}}') (r_{{{nnn}}} - r_{{{nnn}}}') }}{{ || r - r'||^5 }}$"
        )
        ax_complete.set_title(rf"$\Pi(r, r')_{{{mmm},{nnn}}}$")



        for ax in [ax_I, ax_cross, ax_complete]:
            ax.set_xlabel('X')
            ax.set_xlim([-1,0])
            ax.set_ylim([0,1])
            ax.set_zlim([0,1])
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.view_init(elev=42, azim=129, roll=0)

        fig.set_tight_layout(True)
"""



# Plot the dot product with the
vec = np.zeros((3, n, n))
for i in range(n):
    for j in range(n):
        vec[:,i,j] = np.dot(TotalPiMatrix[:,:,i,j], stnnorm)


fig = plt.figure()
axn1 = fig.add_subplot(131, projection='3d')
axn2 = fig.add_subplot(132, projection='3d')
axn3 = fig.add_subplot(133, projection='3d')
axns = [axn1, axn2, axn3]

for idir in range(3):
    scal = vec[idir, :, :]

    # Get the positive values
    posscal = np.NaN*scal
    negscal = np.NaN*scal
    posscal[scal > 0] = scal[scal > 0]
    negscal[scal < 0] = scal[scal < 0]

    logpos = np.log10(posscal)
    logneg = np.log10(-negscal)

    scale = np.max([np.abs(np.nanmax(logpos)),
                    np.abs(np.nanmin(logpos)),
                    np.abs(np.nanmax(logneg)),
                    np.abs(np.nanmin(logneg))
                    ])


    normlogpos = (1+(logpos/scale))/2
    normlogneg = (1+(logneg/scale))/2


    vrange = np.arange(scale-8, scale + 1)
    nints = len(vrange)

    # Create custom colorbar
    ni = nints - 1  # Number of discrete colors you want
    oranges = plt.cm.Oranges(np.linspace(0, 1, ni))  # Sample n colors from a continuous colormap
    blues = plt.cm.Blues(np.linspace(0, 1, ni))  # Sample n colors from a continuous colormap
    discrete_orange = LinearSegmentedColormap.from_list("orange_disc", oranges, N=ni)
    discrete_blue = LinearSegmentedColormap.from_list("orange_disc", blues, N=ni)

    axns[idir].plot_surface(x, y, z, facecolors=discrete_orange(normlogpos))
    axns[idir].plot_surface(x, y, z, facecolors=discrete_blue(normlogneg))




    ax = axns[idir]
    ax.set_xlabel('X')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev=42, azim=129, roll=0)

plt.show()