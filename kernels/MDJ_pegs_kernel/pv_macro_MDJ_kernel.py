# Macro for paraview to visualise the PEGS kernel
import os
from paraview.simple import *
import numpy as np

# Setup some directories
mdir = "/Users/eaton/Documents/Princeton/PEGS_kernels/"
ddir = f"{mdir}/kernels/MDJ_kernel/ensight/"

# Load in the wanted procs
procs = np.loadtxt(f'{ddir}/wanted_procs.txt').astype(int)
data     = {}       # for querying the object by its proc id
allprocs = []       # for all the procs
for iproc in procs:
    d = EnSightReader(registrationName=f'reg1_proc{iproc}.case', CaseFileName=f'{ddir}/reg1_proc{iproc}.case')
    data[iproc] = d
    allprocs.append(d)

default_render_view = GetActiveViewOrCreate('RenderView')

# Group all the procs and display
ProcsGrouped = GroupDatasets(registrationName='ProcsGrouped', Input=allprocs)
ProcsGrouped_Display = Show(ProcsGrouped, default_render_view,  'UnstructuredGridRepresentation')
ProcsGrouped_Display.Representation = 'Surface'


# Add the coastlines:
Coastlines = XMLPolyDataReader(registrationName='Coastlines',
                               FileName=[f'{mdir}/Figures/data_for_plots/World_Coastlines_and_Lakes.vtp'])
# Clip to the Tohoku region
clip = Clip(registrationName='Clip', Input=Coastlines)
clip.Invert = 0
clip.ClipType.Normal = [0.0, 0.0, 1.0]

clip = Clip(registrationName='Clip', Input=clip)
clip.ClipType.Origin = [0.08123012750505833, 0.004245012998580933, 0.5002922564744949]
clip.ClipType.Normal = [1.0, 0.0, 0.0]

clip = Clip(registrationName='Clipped_coastlines', Input=clip)
clip.Invert = 0
clip.ClipType.Normal = [0.0, 1.0, 0.0]

Coastlines_clipped_display = Show(clip, default_render_view, 'UnstructuredGridRepresentation')

# Hide the main surface
Hide(ProcsGrouped, default_render_view)

# Lets start thresholding:
#maxamp = 1e-14
#minamp = 1e-23
#scalebar = 1e-22

maxamp = 1e-20
minamp = 1e-23
scalebar = 1e-22

knl_name = 'rhonotprime'

# Positive threshold
"""thr_pos = Threshold(registrationName=f'thr_{knl_name}_pos', Input=ProcsGrouped)
thr_pos.Scalars = ['POINTS', 'rho_kernel']
thr_pos.LowerThreshold = minamp
thr_pos.UpperThreshold = maxamp

# Negative threshold:
thr_neg = Threshold(registrationName=f'thr_{knl_name}_neg', Input=ProcsGrouped)
thr_neg.Scalars = ['POINTS', 'rho_kernel']
thr_neg.LowerThreshold = -maxamp
thr_neg.UpperThreshold = -minamp

for thr in [thr_neg, thr_pos]:
    display = Show(thr, default_render_view, 'UnstructuredGridRepresentation')
    ColorBy(display, ('POINTS', knl_name))
    Hide(thr, default_render_view)
"""


clip_pos = Clip(registrationName=f'{knl_name}_pos', Input=ProcsGrouped)
clip_pos.ClipType = 'Scalar'
clip_pos.Scalars = ['POINTS', knl_name]
clip_pos.Value = -minamp

clip_neg = Clip(registrationName=f'{knl_name}_neg', Input=ProcsGrouped)
clip_neg.ClipType = 'Scalar'
clip_neg.Scalars = ['POINTS', knl_name]
clip_neg.Value = minamp
clip_neg.Invert = 0


Show(clip_pos, default_render_view)
Show(clip_neg, default_render_view)

# ---- NEXT WE CREATE A SOURCE RECEIVER SLICE ----

def lat_lon_to_xyz(r, lat, lon):
    if lon < 0:
        lon = 360 + lon
    theta = np.deg2rad(lon)
    phi = np.pi / 2 - np.deg2rad(lat)
    xcart = r * np.cos(theta) * np.sin(phi)
    ycart = r * np.sin(theta) * np.sin(phi)
    zcart = r * np.cos(phi)
    return xcart, ycart, zcart

# Compute normals for source-receiver locations in paraview:
src_lat = 37.5200
src_lon = 143.0500

# Station MDJ coordinates
stn_lat = 44.61
stn_lon = 129.5908
r = 1

# Vectors to src and receiver
n_src = lat_lon_to_xyz(r, src_lat, src_lon)
n_stn = lat_lon_to_xyz(r, stn_lat, stn_lon)

# Compute source receiver plane
plane_normal = np.cross(n_src, n_stn)
plane_origin = (np.array(n_stn) + np.array(n_src))/2

# Slice the model in this plane
src_rec_slice = Slice(registrationName=f"src_rec_slice", Input=ProcsGrouped)
src_rec_slice.SliceType.Origin = plane_origin
src_rec_slice.SliceType.Normal = plane_normal

# we want to only have the part of the src-rec slice that is close
# to the src and receiver - 15 % buffer around src-rec arc length
buffsize = 0.15
d = np.array(n_stn) - np.array(n_src)
L = np.linalg.norm(d)
d_hat = d / L

# Create the two vectors that defined the buffer points
n_buffer1 = np.array(n_src) - buffsize * L * d_hat
nplanestn1 =  np.cross(plane_normal, n_buffer1)

n_buffer2 = np.array(n_stn) + buffsize * L * d_hat
nplanestn2 =  np.cross(plane_normal, n_buffer2)

# Clip the plane on each side:
# -- buffer 1
src_rec_slicebuf1 = Clip(registrationName=f"src_rec_slicebuf1", Input=src_rec_slice)
src_rec_slicebuf1.ClipType  = 'Plane'
src_rec_slicebuf1.ClipType.Origin = [0,0,0]
src_rec_slicebuf1.ClipType.Normal = nplanestn1
src_rec_slicebuf1.Invert = 0
# -- buffer 2
src_rec_slicebuf2 = Clip(registrationName=f"src_rec_slicebuf2", Input=src_rec_slicebuf1)
src_rec_slicebuf2.ClipType  = 'Plane'
src_rec_slicebuf2.ClipType.Origin = [0,0,0]
src_rec_slicebuf2.ClipType.Normal = nplanestn2

# Probably want to cut down this slice to top ~ 250 km
topdepth = 250
src_rec_sliceUPPER = Clip(registrationName=f"src_rec_slice_{topdepth}", Input=src_rec_slicebuf2)
src_rec_sliceUPPER.ClipType  = 'Sphere'
src_rec_sliceUPPER.ClipType.Center = [0,0,0]
src_rec_sliceUPPER.ClipType.Radius = (6371-topdepth)/6371
src_rec_sliceUPPER.Invert = 0


# Change the display colour of the slice:
src_rec_slice_Display = GetDisplayProperties(src_rec_sliceUPPER, view=default_render_view)
ColorBy(src_rec_slice_Display, ('POINTS', knl_name))
src_rec_slice_Display.SetScalarBarVisibility(default_render_view, True)

# Get look up table for rho_kernel and rescale to these bounds
rho_kernelLUT = GetColorTransferFunction(knl_name)
rho_kernelLUT.RescaleTransferFunction(-scalebar, scalebar)

# Apply the colorbar name:
rho_kernelLUT.ApplyPreset('broc', True)

# get opacity transfer function/opacity map for 'rho_kernel'
# and rescale transfer function
rho_kernelPWF = GetOpacityTransferFunction(knl_name)
rho_kernelPWF.RescaleTransferFunction(-scalebar, scalebar)
rho_kernelTF2D = GetTransferFunction2D(knl_name)






# # Try to orient the camera:
#
# # ======== Define your three points in the plane ========
# # Replace these with your actual coordinates
# source        = np.array(n_src)
# receiver      = np.array(n_stn)
# surface_point = np.array(n_buffer1)
#
# # ======== Compute vectors and normal ========
# v1 = receiver - source
# v2 = surface_point - source
# normal = np.cross(v1, v2)
# normal /= np.linalg.norm(normal)
#
# # ======== Set camera ========
# focal_point = source  # or use the centroid of the plane
# camera_distance = 1.0  # adjust for zoom level
# camera_position = focal_point + camera_distance * normal
#
# # view-up direction in the plane (ensure it's not aligned with the normal)
# view_up = v1 / np.linalg.norm(v1)
#
# # ======== Set ParaView camera ========
# view = GetActiveViewOrCreate('RenderView')
# view.CameraPosition = camera_position.tolist()
# view.CameraFocalPoint = focal_point.tolist()
# view.CameraViewUp = view_up.tolist()
#
# # Render the view
# Render()





# perpendicular vector to src ppint:
srcpoint = np.array(n_src) * 6351/6371
north = np.cross(n_src,plane_normal)


p1 = srcpoint + north/50
p2 = srcpoint - north/50

vec = p2 - p1

import numpy as np
from scipy.spatial.transform import Rotation as R



# Plane normal (must be a unit vector)
n = plane_normal / np.linalg.norm(plane_normal)

# Define rotation: 10 degrees about the plane normal
angle_deg = 10
rotation = R.from_rotvec(np.deg2rad(angle_deg) * n)
# Apply rotation
v_rotated = rotation.apply(vec)
print(v_rotated)

p3 = srcpoint + v_rotated*2
p4 = srcpoint - v_rotated*2




# Define rotation: 10 degrees about the plane normal
angle_deg = -80
rotation = R.from_rotvec(np.deg2rad(angle_deg) * n)
# Apply rotation
v_rotated = rotation.apply(vec)
print(v_rotated)

p5 = srcpoint + v_rotated*2
p6 = srcpoint - v_rotated*2


