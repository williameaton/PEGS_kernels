import os
import numpy as np
from paraview.simple import *
# GENERATE PROC LIST FROM DIRECTROY:
import numpy as np
strain_kernels = True
save_vtkjs = False 
save_vtm   = False
add_slice_src_rec = False
# Setup some directories
wp_dir = "/Users/eaton/Documents/Princeton/PEGS_kernels/Figures/data_for_plots/MDJ_wanted_procs.txt"
if strain_kernels:
    mdir  = "/Users/eaton/Documents/Princeton/PEGS_kernels/strain/MDJ_kernel/quick_look"
    tag_dir = 'PGSS'
else:
    mdir   = "/Users/eaton/Documents/Princeton/PEGS_kernels/kernels/MDJ_kernel/ensight"
    tag_dir = 'PEGS'


# Load in the wanted procs
procs = np.loadtxt(wp_dir).astype(int)
data     = {}       # for querying the object by its proc id
allprocs = []       # for all the procs
for iproc in procs:
    d = EnSightReader(registrationName=f'reg1_proc{iproc}.case', CaseFileName=f'{mdir}/reg1_proc{iproc}.case')
    data[iproc] = d
    allprocs.append(d)


default_render_view = GetActiveViewOrCreate('RenderView')

# Group all the procs and display
ProcsGrouped = GroupDatasets(registrationName='ProcsGrouped', Input=allprocs)
ProcsGrouped_Display = Show(ProcsGrouped, default_render_view,  'UnstructuredGridRepresentation')
ProcsGrouped_Display.Representation = 'Surface'


if add_slice_src_rec:
    def lat_lon_to_xyz(r, lat, lon):
        if lon < 0:
            lon = 360 + lon
        theta = np.deg2rad(lon)
        phi = np.pi / 2 - np.deg2rad(lat)
        xcart = r * np.cos(theta) * np.sin(phi)
        ycart = r * np.sin(theta) * np.sin(phi)
        zcart = r * np.cos(phi)
        return xcart, ycart, zcart
    #
    # Compute normals for source-receiver locations in paraview:
    src_lat = 37.5200
    src_lon = 143.0500
    #
    stn_lat = 44.61
    stn_lon = 129.5908
    r = 1
    #
    n_src = lat_lon_to_xyz(r, src_lat, src_lon)
    n_stn = lat_lon_to_xyz(r, stn_lat, stn_lon)
    #
    # Compute normal
    plane_normal = np.cross(n_src, n_stn)
    plane_origin = (np.array(n_stn) + np.array(n_src))/2
    #
    # Create slice
    src_rec_slice = Slice(registrationName=f"src_rec_slice", Input=ProcsGrouped)
    src_rec_slice.SliceType.Origin = plane_origin
    src_rec_slice.SliceType.Normal = plane_normal






load_coastlines = True
save_coastlines = True

if load_coastlines or save_coastlines:
    # Add the coastlines clipped around Japan:
    Coastlines = XMLPolyDataReader(registrationName='Coastlines',
                                   FileName=[f'/Users/eaton/Documents/Princeton/PEGS_kernels/Figures/data_for_plots/World_Coastlines_and_Lakes.vtp'])
    # Generate first clip:
    clip = Clip(registrationName='Clip', Input=Coastlines)
    clip.Invert = 0
    clip.ClipType.Normal = [0.0, 0.0, 1.0]
    clip.ClipType.Origin = [0.0, 0.0, 0.0]
    #
    origins = [[0.08123012750505833, 0.004245012998580933, 0.5002922564744949],
               [0,0,0],
               [-0.4026201636517875, 0.5933291302822371, 0.3923995908434186],
               [-0.43010589059946086, 0.3121248113813596, 0.5149118853983539],
               [-0.39925526571325587, 0.4296072646975517, 0.4274214757606387],
               [-0.6097499977816229, 0.43052062782319456, 0.4548765029759808],
               [-0.5097907374707544, 0.664445959524419, 0.6342947330188587]]
    #
    normals = [[1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.3520777091514397, 0.6417745501495115, -0.681297815568982],
               [-0.4040458873841895, 0.9036852977185497, -0.14177377604923735],
               [1,0,0],
               [0.9526646395492466, 0.010108531751109442, 0.30385506765288767],
               [0.033064190033996764, 0.9993401857938649, -0.015031712972228956]]
    #
    inverts = [1, 0, 1, 0, 1, 0, 1]
    #
    nclips = 6
    #
    for iclip in range(nclips):
        clip = Clip(registrationName='Clip', Input=clip)
        origin = origins[iclip]
        if origin == origin: # non a nan
            clip.ClipType.Origin = origin
        clip.ClipType.Normal = normals[iclip]
        clip.Invert = inverts[iclip]
    Coastlines_clipped_display = Show(clip, default_render_view,'UnstructuredGridRepresentation')
    #
    if save_coastlines:
        SaveData(f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/coastline.pvd', proxy=clip, CellDataArrays=['plates'])


# Hide the main surface
Hide(ProcsGrouped, default_render_view)

# Thresholding:
if strain_kernels:
    posthresh = 1e-3
else:
    posthresh =  1e-23

negthresh    = -posthresh
cb_max_scale = posthresh * 10
# This loop loads each of the kernels one at a time (cropped by the threshold)
# it saves the to VTM files which can then be read back in to export as a vtkjs
knls = ['alpha', 'beta', 'rho', 'kappa', 'mu', 'rhonotprime']
for knl_name in knls:
    #
    # Positive Clip:
    PosClip = Clip(registrationName=f'{knl_name}_POS', Input=ProcsGrouped)
    PosClip.ClipType = 'Scalar'
    PosClip.Scalars = ['POINTS', knl_name]
    PosClip.Value = posthresh
    PosClip.Invert = 0
    # NegativeClip:
    NegClip = Clip(registrationName=f'{knl_name}_NEG', Input=ProcsGrouped)
    NegClip.ClipType = 'Scalar'
    NegClip.Scalars = ['POINTS', knl_name]
    NegClip.Value = negthresh
    NegClip.Invert = 1
    # RE-color:
    LUT = GetColorTransferFunction(knl_name)
    LUT.RescaleTransferFunction(-cb_max_scale, cb_max_scale)
    LUT.ApplyPreset('broc', True)
    PWF = GetOpacityTransferFunction(knl_name)
    PWF.RescaleTransferFunction(-cb_max_scale, cb_max_scale)
    PWF.Points = [-cb_max_scale, 1.0, 0.5, 0.0, -0.00010000014282375569, 1.0, 0.5, 0.0, 2.8985426979488693e-05, 1.0, 0.5, 0.0, 9.999986000443911e-05, 1.0, 0.5, 0.0, cb_max_scale, 1.0, 0.5, 0.0]
    TF2D = GetTransferFunction2D(knl_name)
    TF2D.RescaleTransferFunction(-cb_max_scale, cb_max_scale, 0.0, 1.0)
    #
    Show(PosClip, default_render_view)
    Show(NegClip, default_render_view)
    #
    # Save the two clips:
    if save_vtm: 
        SaveData(f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/{knl_name}_POS.vtm', proxy=PosClip, ChooseArraysToWrite=1,PointDataArrays=[knl_name],Assembly='Hierarchy')
        SaveData(f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/{knl_name}_NEG.vtm', proxy=NegClip, ChooseArraysToWrite=1,PointDataArrays=[knl_name],Assembly='Hierarchy')
    #
    # Delete the clips
    Delete(PosClip)
    del PosClip
    Delete(NegClip)
    del NegClip

# Delete the procs grouped:
Delete(ProcsGrouped)
del ProcsGrouped
for iproc in procs:
    d = FindSource(f'reg1_proc{iproc}.case')
    Delete(d)
    del d


for knl_name in knls:
    # Read in the pos and neg bits:
    NEGvtm = XMLMultiBlockDataReader(registrationName=f'{knl_name}_NEG.vtm', FileName=[f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/{knl_name}_NEG.vtm'])
    POSvtm = XMLMultiBlockDataReader(registrationName=f'{knl_name}_POS.vtm', FileName=[f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/{knl_name}_POS.vtm'])
    # Read coastline
    coastlinepvd = PVDReader(registrationName='coastline.pvd', FileName=f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/coastline.pvd')
    #
    NEGvtmDisplay = Show(NEGvtm, default_render_view, 'UnstructuredGridRepresentation')
    POSvtmDisplay = Show(POSvtm, default_render_view, 'UnstructuredGridRepresentation')
    CoastDisplay  = Show(coastlinepvd, default_render_view)
    #
    LUT = GetColorTransferFunction(knl_name)
    LUT.RescaleTransferFunction(-cb_max_scale, cb_max_scale)
    LUT.ApplyPreset('broc', True)
    PWF = GetOpacityTransferFunction(knl_name)
    PWF.RescaleTransferFunction(-cb_max_scale, cb_max_scale)
    PWF.Points = [-cb_max_scale, 1.0, 0.5, 0.0, -0.00010000014282375569, 1.0, 0.5, 0.0, 2.8985426979488693e-05, 1.0, 0.5, 0.0, 9.999986000443911e-05, 1.0, 0.5, 0.0, cb_max_scale, 1.0, 0.5, 0.0]
    TF2D = GetTransferFunction2D(knl_name)
    TF2D.RescaleTransferFunction(-cb_max_scale, cb_max_scale, 0.0, 1.0)
    # Final thing is to trick the coastlines colorbar
    LUTcoastlines = GetColorTransferFunction('plates')
    LUTcoastlines.ApplyPreset('AllBlack', True)
    #
    # Export to VTKJS:
    if save_vtkjs:
        ExportView(f'/Users/eaton/Documents/Personal/web/web/glance/{tag_dir}/{knl_name}.vtkjs', view=default_render_view, ParaViewGlanceHTML='')
    #
    #


    Delete(coastlinepvd)
    del coastlinepvd
    Delete(NEGvtm)
    del NEGvtm
    Delete(POSvtm)
    del POSvtm





