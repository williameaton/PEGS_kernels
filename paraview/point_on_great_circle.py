# Re-importing necessary libraries
import numpy as np

# Redefine constants and utility functions after the reset
deg_to_rad = np.pi / 180
rad_to_deg = 180 / np.pi

# Redefine spherical to Cartesian and Cartesian to spherical conversions
def spherical_to_cartesian(theta, phi):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def cartesian_to_spherical(cartesian):
    x, y, z = cartesian
    theta = np.arccos(z)
    phi = np.arctan2(y, x)
    return theta, phi

# Redefine rotation matrix computation
def rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    ux, uy, uz = axis
    return np.array([
        [cos_a + ux**2 * (1 - cos_a), ux*uy*(1 - cos_a) - uz*sin_a, ux*uz*(1 - cos_a) + uy*sin_a],
        [uy*ux*(1 - cos_a) + uz*sin_a, cos_a + uy**2 * (1 - cos_a), uy*uz*(1 - cos_a) - ux*sin_a],
        [uz*ux*(1 - cos_a) - uy*sin_a, uz*uy*(1 - cos_a) + ux*sin_a, cos_a + uz**2 * (1 - cos_a)]
    ])

def get_plane(epidist):
    # Reinitialize given spherical coordinates (latitude, longitude)
    lat1, lon1 = 37.5200, 143.0500
    lat2, lon2 = 44.6100, 129.5908
    # Convert to colatitude (theta) and longitude (phi) in radians
    theta1 = (90 - lat1) * deg_to_rad
    phi1 = lon1 * deg_to_rad
    theta2 = (90 - lat2) * deg_to_rad
    phi2 = lon2 * deg_to_rad
    # Convert points to Cartesian coordinates
    P1 = spherical_to_cartesian(theta1, phi1)
    P2 = spherical_to_cartesian(theta2, phi2)
    # Compute the normal vector of the great circle
    N = np.cross(P1, P2)
    N /= np.linalg.norm(N)  # Normalize the normal vector
    # Epicentral distance of 12 degrees in radians
    delta_12 = epidist * deg_to_rad
    # Compute the rotation matrix for 12 degrees
    R_12 = rotation_matrix(N, delta_12)
    P_rotated_12 = np.dot(R_12, P1)  # Rotated point for 12 degrees in Cartesian coordinates
    # Convert back to spherical coordinates
    theta_rot_12, phi_rot_12 = cartesian_to_spherical(P_rotated_12)
    # Convert back to latitude and longitude
    lat_rot_12 = 90 - theta_rot_12 * rad_to_deg
    lon_rot_12 = phi_rot_12 * rad_to_deg
    # Ensure longitude is in the range [-180, 180]
    if lon_rot_12 > 180:
        lon_rot_12 -= 360
    elif lon_rot_12 < -180:
        lon_rot_12 += 360
    # Convert the spherical coordinates for 12 degrees to Cartesian
    x_rot_12, y_rot_12, z_rot_12 = spherical_to_cartesian((90 - lat_rot_12) * deg_to_rad, lon_rot_12 * deg_to_rad)
    return np.array([x_rot_12, y_rot_12, z_rot_12])

def lat_lon_to_xyz(r, lat, lon):
    if lon < 0:
        lon = 360 + lon
    theta = np.deg2rad(lon)
    phi = np.pi / 2 - np.deg2rad(lat)
    xcart = r * np.cos(theta) * np.sin(phi)
    ycart = r * np.sin(theta) * np.sin(phi)
    zcart = r * np.cos(phi)
    return xcart, ycart, zcart



# Create the source-receiver great circle slice
# Compute normals for source-receiver locations in paraview:
src_lat = 37.5200
src_lon = 143.0500

stn_lat = 44.61
stn_lon = 129.5908

r = 1

p_src = lat_lon_to_xyz(r, src_lat, src_lon)
p_stn = lat_lon_to_xyz(r, stn_lat, stn_lon)

# Compute normal
plane_normal = np.cross(p_src, p_stn)
plane_origin = (np.array(p_stn) + np.array(p_src))/2

# Create slice
src_rec_slice = Slice(registrationName=f"src_rec_slice", Input=ProcsGrouped)
src_rec_slice.SliceType.Origin = plane_origin
src_rec_slice.SliceType.Normal = plane_normal


# Cut each side:
leftclip = Clip(registrationName='LeftClip', Input=src_rec_slice)
p1 = get_plane(15)
plane_side1 = np.cross(plane_normal, p1)
leftclip.ClipType.Origin = p1
leftclip.ClipType.Normal = plane_side1

rightclip = Clip(registrationName='RightClip', Input=leftclip)
p2 = get_plane(357)
plane_side2 = np.cross(plane_normal, p2)
rightclip.ClipType.Origin = p2
rightclip.ClipType.Normal = plane_side2
rightclip.Invert = 0

# Clip the bottom out
sphereclip = Clip(registrationName='SphereClip', Input=rightclip)
sphereclip.ClipType = 'Sphere'
sphereclip.ClipType.Center = [0.0, 0.0, 0.0]
sphereclip.ClipType.Radius = (6371-250)/6371
sphereclip.Invert = 0


# Since the edges are at epicentral distnaces of 15 and -3 we centre it at 6
centre_up = get_plane(6)

focalpoint = p1 + (p2-p1)/2


renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.CameraPosition   = [0.4645059120918507, 2.7374856221850714, -3.0016418626985515]
renderView1.CameraViewUp     = centre_up
renderView1.CameraFocalPoint = focalpoint
renderView1.CameraViewAngle = 2.5
