import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import cartopy.crs as crs
import numpy as np

# Plots a map of stations based on a SPECFEM STATIONS file
def plot_stn_map_from_file(fname, extent=None, add_ocean=False, add_land=False,
                           add_coastlines=False, coastline_resolution='10m', stnmarkersize=40,
                           labels=None, marker_colour='gold'):
    fig = plt.figure(figsize=(6, 4))
    fig.set_tight_layout(True)
    ax = plt.axes(projection=crs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, alpha=0.2)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Load file data:
    stnfile = open(fname, 'r')
    lines = stnfile.readlines()
    stnname = []; stnnet = []; stnlat = []; stnlon = [];
    stnfeat = [stnname, stnnet, stnlat, stnlon]
    nstns = 0
    for l in lines:
        ll = l.split()
        for il in range(4):
            if il>1:
                stnfeat[il].append(float(ll[il]))
            else:
                stnfeat[il].append(ll[il])
        nstns+=1



    # Set extent:
    if extent!=None:
        map_extent = extent
    else:
        buffer = 2
        min_lat = min(stnlat)
        max_lat = max(stnlat)
        min_lon = min(stnlon)
        max_lon = max(stnlon)
        map_extent = [min_lon - buffer, max_lon + buffer, min_lat-buffer, max_lat+buffer]
    ax.set_extent(map_extent)

    if add_ocean:
        ax.add_feature(cfeature.OCEAN, alpha=0.3)
    if add_land:
        ax.add_feature(cfeature.LAND, alpha=0.3)
    if add_coastlines:
        ax.coastlines(coastline_resolution, alpha=0.7)

    # Plot stations
    ax.scatter(x=stnlon, y=stnlat, c=marker_colour, edgecolors='k', s=stnmarkersize, transform=crs.PlateCarree(),
                     marker='^', zorder=3)


    # Add labels:
    if labels!=None:
        if len(labels)==3:
            stnname = labels[0]
            labx    = labels[1]
            laby    = labels[2]
        elif len(labels)==2:
            labx    = labels[0]
            laby    = labels[1]
        else:
            raise ValueError('Format for labels should be either a list of len 2 or 3: \n  -- len 2: [[lon buffers],[lat buffers]], \n  -- len 3: [[names], [lon buffers],[lat buffers]]')

    else:
        labx = 0.1 + np.zeros(nstns)
        laby = 0.1 + np.zeros(nstns)

    for istn in range(nstns):
        ax.text(x=stnlon[istn]+labx[istn], y=stnlat[istn]+laby[istn], s=stnname[istn], transform=crs.PlateCarree())


    return fig, ax





