# Testing the rotation of QSSP acceleration data to produce strain data as required
import math
from obspy.geodetics import gps2dist_azimuth
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(f"../../../classes")
from finite_difference import FD, FdStation
from wetools.plotting import setup_we_mpl, subplots_hide_xaxes




def rotate_strain(hxyz, baz):
    # Taken from Zhang supp data.
    ''' rotate hxyz to hrtz according to back-azimuth
    hxyz: hxx, hyy, hxy, hyx, hxz, hyz
    hrtz: hrr, htt, hrt, htr, hrz, htz
    '''
    alpha = 270-baz    # in degree
    ca = math.cos(alpha/180*math.pi)
    sa = math.sin(alpha/180*math.pi)
    hrr = ca*ca*hxyz[:,0]  + sa*sa*hxyz[:,1] + 2*ca*sa*hxyz[:,2]
    htt = sa*sa*hxyz[:,0]  + ca*ca*hxyz[:,1] - 2*ca*sa*hxyz[:,2]
    hrt = -ca*sa*hxyz[:,0] + ca*sa*hxyz[:,1] + (ca*ca-sa*sa)*hxyz[:,2]
    htr = hrt
    hrz = ca*hxyz[:,4]  + sa*hxyz[:,5]
    htz = -sa*hxyz[:,4] + ca*hxyz[:,5]
    hrtz = np.stack((hrr,htt,hrt,htr,hrz,htz), axis=1)
    return hrtz

chls      = ['ZZ', 'NN', 'EE', 'ZN', 'EN', 'ZE']
pltchls   = ['++', 'ZZ', 'XX', 'XX', 'RZ', 'TZ']
nplotchls = len(pltchls)
dx = 5000

zhang_ddir = '../../data/ZHANG_data/'

# Setup
setup_we_mpl()
fd = FD(from_readme=True, readme_fname='qssp_zhang_stns_plus_PEGSstns/README')

fig, ax = plt.subplots(6,6, figsize=(20,10.5))
#fig.set_tight_layout(True)


# Compute strain from centred FD3 of QSSP acc
# load the acceleration data
fd.load_qssp_grav_acc(fpath=f'{zhang_ddir}/Data_Zhang2023/HalfSpace/')
fd.strain_from_acc(method=f'centred_FD3', dx=dx, code='qssp')


stnctr = 0
for istn in fd.stations:
    qs   = istn.qssp

    time = istn.qssp["time"]
    fs = np.mean(time[1:]-time[:-1])

    for i in ['EE', 'NN', 'EN', 'NE', 'ZE', 'ZN']:
        dd = qs[f'{i}.STRAIN']
        dd = cumtrapz(dd, dx=1.0/fs, initial=0)
        dd = cumtrapz(dd, dx=1.0/fs, initial=0)
        qs[f'{i}.STRAIN'] = dd


    # This now needs to be rotated to the orientation used in Zhang paper:
    h_ENZ = np.stack((qs['EE.STRAIN'],
                      qs['NN.STRAIN'],
                      qs['EN.STRAIN'],
                      qs['NE.STRAIN'],
                      qs['ZE.STRAIN'],
                      qs['ZN.STRAIN']), axis=1)
    _, faz, baz = gps2dist_azimuth(37.5200, 143.0500, istn.clat, istn.clon,
                                   a=6371e3, f=0.0)

    h_RTZ = rotate_strain(h_ENZ,baz)
    # h++ = (hrr-htt)/2
    h_plus = 0.5*(h_RTZ[:,0]-h_RTZ[:,1])
    # hrr+htt+hzz=0
    h_zz = -1.0*(h_RTZ[:,0]+h_RTZ[:,1])
    # h_RTZ: h++,hzz,hxx,hxx,hrz,htz
    h_RTZ[:,0] = h_plus
    h_RTZ[:,1] = h_zz
    h_RTZ*=1e12

    for i in range(6):
        ax[stnctr,i].plot(time, h_RTZ[:,i])



    stnctr+=1


# directly load the strain from files and plot
fd.load_qssp_strain(fpath=f'{zhang_ddir}/Data_Zhang2023/HalfSpace/')

# Plot the QSSP directly supplied from paper:
stnctr = 0
for istn in fd.stations:
    for ichl in range(nplotchls):
        q = istn.qssp
        ax[stnctr,ichl].plot(q["strain_time"], q[istn.name + f'_{pltchls[ichl]}'], color='r')

        if stnctr==0:
            ax[stnctr, ichl].set_title(f"Channel: {pltchls[ichl]}", weight='bold')
    stnctr += 1



ax[-1, 2].legend(['Finite difference of RCVR*_grav_*.dat', 'RCVR*_grav.h.dat'],
                 ncol=2,
                 loc='lower center',
                 bbox_to_anchor=(0.5, -0.5))

for i in range(6):
    ax[i, 0].set_ylabel(f"Station {i+1}",weight='bold')
    ax[i, 0].yaxis.set_label_coords(-0.2, 0.5)




#subplots_hide_xaxes(ax=ax, nrows=6, ncols=6, keep=[2])
subplots_hide_xaxes(ax=ax, nrows=6, ncols=6)


plt.savefig('compare_QSSP_strain_FD.pdf', format='pdf')