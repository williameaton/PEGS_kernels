import numpy as np
#

def slerp(p0, p1, t):
    """
    Spherical linear interpolation between p0 and p1 by fraction t.
    Assumes p0 and p1 are unit vectors.
    """
    p0 = p0 / np.linalg.norm(p0)
    p1 = p1 / np.linalg.norm(p1)
    omega = np.arccos(np.clip(np.dot(p0, p1), -1, 1))  # angle between p0 and p1
    if np.isclose(omega, 0):
        return p0  # or (p0 + p1)/2 normalized
    sin_omega = np.sin(omega)
    return (np.sin((1 - t) * omega) * p0 + np.sin(t * omega) * p1) / sin_omega

def xyz_to_lat_lon(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = 90 - np.rad2deg(np.arccos(z / r))
    lon = np.rad2deg(np.arctan2(y, x))
    if lon > 180:
        lon -= 360  # keep in [-180, 180]
    return lat, lon

def lat_lon_to_xyz(r, lat, lon):
    if lon < 0:
        lon = 360 + lon
    theta = np.deg2rad(lon)
    phi = np.pi / 2 - np.deg2rad(lat)
    xcart = r * np.cos(theta) * np.sin(phi)
    ycart = r * np.sin(theta) * np.sin(phi)
    zcart = r * np.cos(phi)
    return xcart, ycart, zcart

# Tohoku src
src_lat = 37.5200
src_lon = 143.0500

# MDJ
stn_lat = 44.61
stn_lon = 129.5908

r = 1

# Normal to the station and source
n_src = lat_lon_to_xyz(r, src_lat, src_lon)
n_stn = lat_lon_to_xyz(r, stn_lat, stn_lon)

print('Source: '  , n_src)
print('Receiver: ', n_stn)

# Compute normal for source-receiver plane
plane_normal = np.cross(n_src, n_stn)
plane_origin = (np.array(n_stn) + np.array(n_src))/2

print(plane_normal)

# Compute great-circle midpoint
midpoint_cart = slerp(n_src, n_stn, 0.5)



print("Cartesian midpoint on sphere:", midpoint_cart)

# Normalize vectors
n_src_unit = n_src / np.linalg.norm(n_src)
n_stn_unit = n_stn / np.linalg.norm(n_stn)

# Compute angle between them
dot_product = np.dot(n_src_unit, n_stn_unit)
angle_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))  # arc length in radians

# Halfway along the arc
half_arc_length = angle_rad / 2

# Optionally convert to degrees or km (assuming Earth's radius ~6371 km)
angle_deg = np.rad2deg(half_arc_length)
arc_length_km = half_arc_length * 6371  # for Earth

print(f"Half arc distance: {half_arc_length:.6f} radians")
print(f"                = {angle_deg:.6f} degrees")
print(f"                = {arc_length_km:.2f} km")


mid_lat, mid_lon = xyz_to_lat_lon(*midpoint_cart)

print(f"Midpoint latitude: {mid_lat:.6f}°")
print(f"Midpoint longitude: {mid_lon:.6f}°")



from obspy.taup import TauPyModel

# Load PREM model
model = TauPyModel(model="prem")

# Set epicentral distance
distance_deg = angle_deg * 2

# Get ray paths for P phase
paths = model.get_ray_paths(distance_in_degree=distance_deg, phase_list=["P"], source_depth_in_km=20)

pp = paths[0].path

ax = paths.plot_rays(plot_type='cartesian')

ax.axhline(30, color='r')

import matplotlib.pyplot as plt
plt.show()

# ppth = []
# for ipt in pp:
#
#     pt = slerp(n_src, n_stn, np.rad2deg(ipt[2])/distance_deg)
#
#     # pt is a vector so now normalise for the depth
#     depth = ipt[3]
#     pt *= (1 - ipt[3]/6371)
#
#     ppth.append(pt)
#
# np.savetxt('path.txt', np.array(ppth))
#
