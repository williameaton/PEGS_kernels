import numpy as np
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth

def create_rotation_matrix_NEZ_to_XYZ(lat, lon):

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    Z = np.array([np.cos(lat)*np.cos(lon),
                  np.cos(lat)*np.sin(lon),
                  np.sin(lat)])

    N = np.array([-np.sin(lat)*np.cos(lon),
                  -np.sin(lat)*np.sin(lon),
                  np.cos(lat)])

    E = np.array([-np.cos(lat)*np.sin(lon),
                  np.cos(lat)*np.cos(lon),
                  0])

    # Rotation matrix:
    R = np.array([N, E, Z]) # rotates from XYZ --> NEZ


    # returns matrix that rotates NEZ --> XYZ with R.T
    return R.T