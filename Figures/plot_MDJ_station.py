import cartopy.crs as crs
import matplotlib.pyplot as plt
from plot_station_map_from_file import plot_stn_map_from_file

from colour_schemes import hex

HighContrast = hex.Hexes().hex['HighContrast_PT']

def plot_MDJ(dpath):
    # Source location
    src_lat = 37.5200
    src_lon = 143.05

    # Plotting values:
    stnmarkersize = 40

    # Generate map using PEGS stations;
    fig, ax = plot_stn_map_from_file(fname=dpath,
                           add_ocean=True,
                           add_land=True,
                           add_coastlines=True,
                           stnmarkersize=stnmarkersize,
                           labels=[[0.5], [0.5]],
                           marker_colour=HighContrast[1],
                           extent=[127.5, 145, 30, 47.5])

    # Plot the source
    ax.scatter(x=src_lon, y=src_lat, c=HighContrast[2], edgecolors='k', s=150,transform=crs.PlateCarree(), marker='*',zorder=4 )
    ax.text(x=src_lon - 1.5, y=src_lat - 1.5, s='Tohoku\n2011', transform=crs.PlateCarree())

    #ax.legend(['PEGS stations','Tohoku 2011 source'], loc='upper left')

    fig.set_tight_layout(True)

    return fig

if __name__ == "__main__":
    path = './data_for_plots/MDJ_COORDINATES'
    fig = plot_MDJ(path)
    plt.savefig("./pdfs/MDJ_station.pdf",format='pdf')

