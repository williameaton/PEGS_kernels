import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz, trapz
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

def rotation_matrix(elat,elon,rlat,rlon, geoco=1):
    # Rotation matrix for RTZ to NEZ
    # Event lat, event lon, receiver lat, receiver lon, flattening

    thetas = np.pi/2  - np.arctan(geoco * np.tan(np.deg2rad(elat)))
    thetar = np.pi/2  - np.arctan(geoco * np.tan(np.deg2rad(rlat)))
    phis = np.deg2rad(elon)
    phir = np.deg2rad(rlon)

    theta = np.arccos(np.cos(thetar)*np.cos(thetas) + np.sin(thetar)*np.sin(thetas)*np.cos(phir - phis))
    sint = np.sin(theta)

    if np.abs(sint)<1e-9:
        sint=1e-9

    rot1 = (np.sin(thetar)*np.cos(thetas) - np.cos(thetar)*np.sin(thetas)*np.cos(phir - phis)) / sint
    rot2 = np.sin(thetas)*np.sin(phir - phis)/sint

    R = np.zeros((3,3))
    R[0,0] = rot1
    R[0,1] = rot2
    R[1,0] = -rot2
    R[1,1] = rot1
    R[2,2] = 1.0
    return R

def rotate_strain(hxyz, baz):
    import math
    ''' rotate hxyz to hrtz according to back-azimuth
    hxyz: hxx, hyy, hxy, hyx, hxz, hyz
    hrtz: hrr, htt, hrt, htr, hrz, htz
    
    WE: this is using the rotation matirix as follows: 
       cos a    sin a   0
      -sin a    cos a   0
        0        0      1 
    '''
    alpha = 270-baz    # in degree
    ca = math.cos(alpha/180*math.pi)
    sa = math.sin(alpha/180*math.pi)
    hrr = ca*ca*hxyz[:,0] + sa*sa*hxyz[:,1] + 2*ca*sa*hxyz[:,2]
    htt = sa*sa*hxyz[:,0] + ca*ca*hxyz[:,1] - 2*ca*sa*hxyz[:,2]
    hrt = -ca*sa*hxyz[:,0] + ca*sa*hxyz[:,1] + (ca*ca-sa*sa)*hxyz[:,2]
    htr = hrt
    hrz = ca*hxyz[:,4] + sa*hxyz[:,5]
    htz = -sa*hxyz[:,4] + ca*hxyz[:,5]

    hrtz = np.stack((hrr,htt,hrt,htr,hrz,htz), axis=1)

    return hrtz



def Rmat(elat,elon,rlat,rlon):
    from obspy.geodetics import gps2dist_azimuth

    _, faz, baz = gps2dist_azimuth(elat, elon, rlat, rlon,
                                   a=6371e3, f=0.0)

    cosb = np.cos(np.radians(baz))
    sinb = np.sin(np.radians(baz))

    R_lhs = np.zeros((3,3))
    R_lhs[0,0] = -cosb
    R_lhs[1,1] = -cosb
    R_lhs[0,1] = -sinb
    R_lhs[1,0] =  sinb
    R_lhs[2,2] = 1.0

    R_rhs = np.zeros((3,3))
    R_rhs[0,0] = -cosb
    R_rhs[1,1] = cosb
    R_rhs[0,1] = -sinb
    R_rhs[1,0] = -sinb
    R_rhs[2,2] = 1.0

    return R_lhs, R_rhs


#obspy_rotate(elat=0, elon=0, rlat=-40, rlon=20)



from reproduce_Zhang_Figs import reproduce_ZhangFig2,reproduce_ZhangFig3, plot_zhang

zhang_ddir    = '../data/ZHANG_data/Data_Zhang2023'
spfm_rootdir  = './archive_specfem/'

fig, ax     = reproduce_ZhangFig2(fpath=zhang_ddir)
figRZ, axRZ = reproduce_ZhangFig3(fpath=zhang_ddir)

figZZ, axZZ = plot_zhang(fpath=zhang_ddir, index=2)
figRT, axRT = plot_zhang(fpath=zhang_ddir, index=3)
figTZ, axTZ = plot_zhang(fpath=zhang_ddir, index=6)



elat =  37.5200
elon =  143.0500



stn_coords =  {"HS1":  [38.5 ,  140.5],
               "HS2":  [37.75,  139.5],
               "HS3":  [37.0 ,  140.5],
               "HS4":  [45.0 ,  135.3],
               "HS5":  [45.0 ,  134.0],
               "HS6":  [43.7 ,  134.0]}

stn_cuts =    {"HS1":  32,
               "HS2":  40,
               "HS3":  30,
               "HS4":  135,
               "HS5":  144,
               "HS6":  132}

# Load the data for the homogenous halfspace:
datadir = './halfspace_data/no_time_int'

colctr = 0
rowctr = 0



for stnno in range(6):
    stn = f'HS{stnno+1}'

    # Get station coordinates
    rlat = stn_coords[stn][0]
    rlon = stn_coords[stn][1]


    _, faz, baz = gps2dist_azimuth(elat, elon, rlat, rlon,
                                   a=6371e3, f=0.0)

    R = np.zeros((3,3))
    theta = 270 - baz
    costh = np.cos(np.deg2rad(theta))
    sinth = np.sin(np.deg2rad(theta))

    R[0,0] =  costh
    R[0,1] =  sinth
    R[1,0] = -sinth
    R[1,1] =  costh
    R[2,2] =  1


    # Load specfem_homogenous:
    clrs = ['crimson', 'teal', 'darkgreen','purple', ]
    ictr = 0

    for method in ['rho2000_FD_perfect_sphere/1_step_raw',
                   'rho2000_FD_perfect_sphere/FD',
                   'rho2000_FD_perfect_sphere/twostep']:

        clr = clrs[ictr]
        ictr +=1

        spfmx_dir = spfm_rootdir + method
        label_spfmx = method

        if 'globe' in spfm_rootdir:
            nn_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HNN.PGRAV.sem.ascii')[:1100,:]
            ee_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HEE.PGRAV.sem.ascii')[:1100,:]
            zz_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HZZ.PGRAV.sem.ascii')[:1100,:]
            nz_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HNZ.PGRAV.sem.ascii')[:1100,:]
            ez_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HEZ.PGRAV.sem.ascii')[:1100,:]
            ne_spec = np.loadtxt(f'{spfmx_dir}/YY.{stn}_33.HNE.PGRAV.sem.ascii')[:1100,:]
        else:
            nn_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXNN.STRAIN.sem.ascii')[:, :]
            ee_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXEE.STRAIN.sem.ascii')[:, :]
            zz_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXZZ.STRAIN.sem.ascii')[:, :]
            nz_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXNZ.STRAIN.sem.ascii')[:, :]
            ez_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXEZ.STRAIN.sem.ascii')[:, :]
            ne_spec = np.loadtxt(f'{spfmx_dir}/{stn}_33.YY.BXNE.STRAIN.sem.ascii')[:, :]


        tspec    = ne_spec[:,0] + 70
        dx       = np.mean(tspec[1:]-tspec[:-1])
        lentime  = len(tspec)


        # ----------------- ROTATION FROM NMSYN METHOD -----------------
        # Build matrix for rotation:
        D = np.zeros((lentime, 3, 3))

        D[:,0,0] = ee_spec[:,1]
        D[:,1,1] = nn_spec[:,1]
        D[:,2,2] = zz_spec[:,1]

        D[:,0,1] = ne_spec[:,1]
        D[:,1,0] = ne_spec[:,1]

        D[:,2,0] = ez_spec[:,1]
        D[:,0,2] = ez_spec[:,1]

        D[:,1,2] = nz_spec[:,1]
        D[:,2,1] = nz_spec[:,1]



        # Rotate backwards from NEZ --> RTZ so R^T D R  instead of R D R^T
        for i in range(lentime):
            D[i,:,:] = R@D[i,:,:]@R.T
        first_int  = cumtrapz(D, dx=dx, initial=0, axis=0)
        second_int = cumtrapz(first_int, dx=dx, initial=0, axis=0)

        second_int *= -1e12 * 1# 5.62/5.31 # zhang uses wrong magnitude


        # Cut time
        tmask = tspec < stn_cuts[stn]
        tspec = tspec[tmask]

        hplus = (second_int[tmask,0,0] - second_int[tmask,1,1])/2

        # Plot ZZ
        axZZ[stnno - rowctr, colctr].plot(tspec, second_int[tmask,2,2], color=clr, label=label_spfmx)

        # Plot RZ
        axRZ[stnno - rowctr, colctr].plot(tspec, second_int[tmask,0,2],  color=clr, label=label_spfmx)

        # Plot TZ
        axTZ[stnno - rowctr, colctr].plot(tspec, second_int[tmask,1,2],  color=clr, label=label_spfmx)

        # Plot RT (cross)
        axRT[stnno - rowctr, colctr].plot(tspec, second_int[tmask,1,0],  color=clr, label=label_spfmx)

        # Plot Hplus
        ax[stnno - rowctr, colctr].plot(tspec, hplus,  color=clr, label=label_spfmx)


        # ----------------- ZHANG ADAPTED METHOD -----------------
        """ee_spec[:,1] = cumtrapz(ee_spec[:,1], dx=dx, initial=0)
        ee_spec[:,1] = cumtrapz(ee_spec[:,1], dx=dx, initial=0)

        zz_spec[:,1] = cumtrapz(zz_spec[:,1], dx=dx, initial=0)
        zz_spec[:,1] = cumtrapz(zz_spec[:,1], dx=dx, initial=0)

        ne_spec[:,1] = cumtrapz(ne_spec[:,1], dx=dx, initial=0)
        ne_spec[:,1] = cumtrapz(ne_spec[:,1], dx=dx, initial=0)

        ez_spec[:,1] = cumtrapz(ez_spec[:,1], dx=dx, initial=0)
        ez_spec[:,1] = cumtrapz(ez_spec[:,1], dx=dx, initial=0)

        nn_spec[:,1] = cumtrapz(nn_spec[:,1], dx=dx, initial=0)
        nn_spec[:,1] = cumtrapz(nn_spec[:,1], dx=dx, initial=0)

        nz_spec[:,1] = cumtrapz(nz_spec[:,1], dx=dx, initial=0)
        nz_spec[:,1] = cumtrapz(nz_spec[:,1], dx=dx, initial=0)

        h_ENZ = 1e12*np.stack((ee_spec[:,1], nn_spec[:,1], ne_spec[:,1], ne_spec[:,1], ez_spec[:,1], nz_spec[:,1]), axis=1)

        h_RTZ = -1*rotate_strain(h_ENZ, baz)
        # h++ = (hrr-htt)/2
        h_plus = 0.5 * (h_RTZ[:, 0] - h_RTZ[:, 1])

        axZZ[stnno - rowctr, colctr].plot(tspec, -1e12*zz_spec[tmask,1], label=label_spfmx, color=clr)

        ax[stnno - rowctr, colctr].plot(tspec, h_plus[tmask], label=label_spfmx, color=clr)
        axRZ[stnno - rowctr, colctr].plot(tspec, h_RTZ[tmask,4], label=label_spfmx, color=clr)

        axRT[stnno - rowctr, colctr].plot(tspec, h_RTZ[tmask,3], label=label_spfmx, color=clr)
        axTZ[stnno - rowctr, colctr].plot(tspec, h_RTZ[tmask,5], label=label_spfmx, color=clr)"""





    # Update column indexes for plot
    if stnno==2:
        colctr = 1
        rowctr = 3


fs = '7'
for i in range(3):
    for j in range(2):
        for axi in [ax, axZZ, axRZ, axTZ, axRT]:
            axi[i,j].legend(fontsize=fs)
            axi[i,j].legend(fontsize=fs)


fig.suptitle('H++')
figZZ.suptitle('HZZ')
figRZ.suptitle('HRZ')
figTZ.suptitle('HTZ')
figRT.suptitle('HRT')

# Set plot limits
for j in range(3):

    axRZ[j, 0].set_xlim([-10,50])
    axRZ[j, 1].set_xlim([-10,150])


# AX RZ:
axRZ[0,0].set_ylim([-1, 4])
axRZ[1,0].set_ylim([-2, 6])
axRZ[2,0].set_ylim([-1, 3])

axRZ[0,1].set_ylim([-4, 6])
axRZ[1,1].set_ylim([-5, 5])
axRZ[2,1].set_ylim([-4, 6])



pdf_name = 'Strain_homogenous_benchmark.pdf'
# Write output to a pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
for fig in [figZZ, fig, figRT, figRZ, figTZ]:
    pdf.savefig(fig)
pdf.close()

