import sys
sys.path.append(f"../../../classes")
from finite_difference import FD, FdStation
import numpy as np
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from wetools.geodetics import dist_lon_at_lat
from obspy.taup.velocity_model import VelocityModel
from obspy.taup.taup_create import build_taup_model
import os

# Create a FD setup:


geoco = 1
earth_rad  = 6371000
source_lat = 37.5200
source_lon = 143.05


fd = FD('qssp_zhang_stns_plus_PEGSstns')

# List of stations we want to use:
# Note that these are based on the Zhang input file to QSSP

names = ['HS1', 'HS2', 'HS3', 'HS4', 'HS5', 'HS6', 'FUK', 'MAJO', 'INU',
         'NE93','BJT', 'MA2', 'XAN', 'ULN', 'SHR', 'MDJ', 'INCN']
nstns = len(names)

lats  = [38.5,     37.75,    37.0,     45.0,     45.0,     43.7 ,
         32.7177, 36.54567, 35.35, 46.927, 40.0183, 59.5756, 34.0313, 47.8651, 44.0563, 44.617, 37.4776 ]
lons  = [140.5,    139.5,    140.5,    135.3,    134.0,    134.0,
         128.7572,  138.2041, 137.029, 122.4302, 116.1679, 150.77, 108.9237, 107.0532, 144.9944, 129.5908, 126.6239 ]

distance_fd = 5000 # metres
timeoffset = 70


# Compute ray theoretical travel times:
filename = "homo_obspy_model.tvel"
vmodel = VelocityModel.read_tvel_file(filename)
taup_model = build_taup_model(filename, output_folder=os.getcwd())
model = TauPyModel(model=f"{filename[:-5]}.npz")



for istn in range(nstns):

    lllat = lats[istn]
    lllon = lons[istn]

    # Compute distance
    dist, faz, baz = gps2dist_azimuth(lat1=source_lat, lon1=source_lon,
                     lat2=lllat, lon2=lllon,
                     a=earth_rad,
                     f=0.0)

    # Convert propagation distance in metres to degrees along great circle:
    prop_deg = 360*dist/(2*np.pi*earth_rad)

    arrivals = model.get_ray_paths(
        source_depth_in_km=20,
        distance_in_degree=prop_deg,
    )

    # P wave arrival time
    cutoff = arrivals[0].time

    # Get the dlat and dlon based on the latitude
    dlon = dist_lon_at_lat(dist=distance_fd, lat=lllat, rad=earth_rad, geoco=geoco)
    dlat = 360*distance_fd/(2*np.pi*earth_rad)

    fd.add_station(FdStation(stnname=names[istn],
                             lat=lllat,
                             lon=lllon,
                             dlat=dlat,
                             dlon=dlon,
                             geoco=geoco,
                             stencil=5,
                             order=2,
                             time_offset=timeoffset,
                             time_cutoff=cutoff
                             ))



fd.write_specfem_stations()


fd.write_readme()