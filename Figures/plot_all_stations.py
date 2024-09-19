import cartopy.crs as crs
import matplotlib.pyplot as plt
from plot_station_map_from_file import plot_stn_map_from_file

from colour_schemes import hex

HighContrast = hex.Hexes().hex['HighContrast_PT']

def plot_all_stations(dpath):
    # Source location
    src_lat = 37.5200
    src_lon = 143.05

    # Strain locations
    strain_stn_lats = [38.5, 37.75, 37.0, 45.0, 45.0, 43.70]
    strain_stn_lon  = [140.5,139.5, 140.5, 135.30, 134.0,134.0 ]


    # Plotting values:
    stnmarkersize = 40

    # Generate map using PEGS stations;
    fig, ax = plot_stn_map_from_file(fname=dpath,
                           add_ocean=True,
                           add_land=True,
                           add_coastlines=True,
                           stnmarkersize=stnmarkersize,
                           labels=[[0.5,-3,0.5,0.5,-4,0.5,0.2,-4, -2],
                                   [0.5,0.5,0.5,0.5,-1.5,0.5,-1.8,1, 1]],
                           marker_colour=HighContrast[1])


    ax.scatter(x=strain_stn_lon, y=strain_stn_lats, c=HighContrast[3],edgecolors='k', s=stnmarkersize,transform=crs.PlateCarree(), marker='^',zorder=3 )

    # Plot the source
    ax.scatter(x=src_lon, y=src_lat, c=HighContrast[2], edgecolors='k', s=150,transform=crs.PlateCarree(), marker='*',zorder=4 )

    # Label strain
    for i in range(6):
        labx=[0.1,-0.8,0.5,0.2,-0.4,-0.2]
        laby=[0.4,0.6,-1.2,0.3,0.6,-1.8]
        ax.text(x=strain_stn_lon[i]+labx[i], y=strain_stn_lats[i]+laby[i], s=str(i+1),transform=crs.PlateCarree())

    ax.legend(['PEGS stations', 'Strain stations', 'Tohoku 2011 source'], loc='upper left')

    fig.set_tight_layout(True)

    return fig

if __name__ == "__main__":
    path = './data_for_plots/PEGS_STATIONS'
    fig = plot_all_stations(path)
    plt.savefig("./pdfs/all_stations.pdf",format='pdf')

