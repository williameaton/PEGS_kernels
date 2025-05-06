import numpy as np
import matplotlib.pyplot as plt

def xyz_to_lat_lon(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = 90 - np.rad2deg(np.arccos(z / r))
    lon = np.rad2deg(np.arctan2(y, x))
    if lon > 180:
        lon -= 360  # keep in [-180, 180]
    return lat, lon

# In this case I am taking the sphere position from placement in paraview
# so that it lines up relatively well with one of the two nodes of the kernel
x = -0.625
y =  0.482
z =  0.605

mid_lat, mid_lon = xyz_to_lat_lon(x, y, z)

print(f"Latitude: {mid_lat:.6f}°")
print(f"Longitude: {mid_lon:.6f}°")

# Radius from centre
rad =  np.sqrt(x**2 + y**2 + z**2)
print('Radius: ', rad*6371)
# using 6335.779


# Sphere radius:
rr = 100        #  km
print(rr/6371)
