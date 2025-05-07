# Plotting the FD time derivatives of the acceleration used in the adjoint source:
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import obspy
from scipy.signal.windows import tukey
from scipy.special import erf
import matplotlib.gridspec as gridspec

import numpy as np
from wetools.funcs import obspy_gen_mpl


def cs(time):
    alpha = 0.15
    window = tukey(len(time), alpha)
    return window

stn = 'ULN'
fpath_master = '../../data/ULN_pegs_kernel'

datadir  = f'{fpath_master}/forward_output/'
adj_path = f'{fpath_master}/adj_source/'
channel = 'Z'

fig = plt.figure(figsize=(4,7.3))

spec = gridspec.GridSpec(ncols=1, nrows=15, figure=fig)

axacc = fig.add_subplot(spec[0:4, 0])
axtap = fig.add_subplot(spec[4:6, 0])
ax2nd = fig.add_subplot(spec[6:10, 0])
axfin = fig.add_subplot(spec[10:14, 0])

axxax = fig.add_subplot(spec[-1, 0])


types =  [f'{channel}acc', f'{channel}grav']

store = [[], []]

for i_trtype in range(len(types)):
    trtype = types[i_trtype]
    fpath = f'./{datadir}/{trtype}/{stn}_proc.sac'

    # Load in the sac:
    trace = obspy.read(fpath)

    # Convert to numpy
    time, acc = obspy_gen_mpl(trace[0])
    dx = np.mean(time[1:] - time[:-1])


    # calculate 5 point finite difference in time
    acc_dx = (-acc[4:] + 8*acc[3:-1] - 8*acc[1:-3] + acc[:-4])/(12*dx)
    acc_2dx =  (-acc[4:] + 16*acc[3:-1] -30*acc[2:-2] + 16*acc[1:-3] - acc[:-4])/(12*dx*dx)

    store[i_trtype].append(acc)
    store[i_trtype].append(acc_dx)
    store[i_trtype].append(acc_2dx)


# REAL DATA:
pegs    = store[1][0] + store[0][0]
d_pegs  = store[1][1] + store[0][1]
dd_pegs = store[1][2] + store[0][2]

# Net values -- acc
axacc.plot(time, pegs, 'thistle')
# 2nd derivative
ax2nd.plot(time[2:-2], dd_pegs, 'thistle')


# Create a best fit polynomial to pegs signal:
deg = 6
p = np.polyfit(x=time, y=pegs, deg=deg)

p_pegs    = time*0
p_pegs_d  = time*0
p_pegs_dd = time*0

for id in range(deg+1):
    p_pegs +=  p[id]*(time**(deg-id))
for id in range(deg):
    p_pegs_d  += (deg-id)*p[id]*(time**(deg-1-id))
for id in range(deg-1):
    p_pegs_dd += (deg-id)*(deg-1-id)*p[id]*(time**(deg-2-id))

# Original quartic fit to acceleration
axacc.plot(time,  p_pegs,    'darkred', linestyle='--')
ax2nd.plot(time,  p_pegs_dd, 'darkred', linestyle='--')

# Create error function
xerf = np.linspace(-4, 7, len(time))
yerf = 0.5 + 0.5 * erf(xerf)
#axt.plot(time, yerf, 'coral')

# Plot inset zoom in of the quartic vs PEGS acceleration
#inset_axes = inset_axes(axacc,
#                        loc = 'center left',
#                        width="50%",
#                        height = "50%",
#                        bbox_to_anchor=(0.03, 0.06, 1, 1),
#                        bbox_transform=axacc.transAxes
#                        )
#inset_axes.tick_params(labelleft=False, labelbottom=False, left=False, bottom=False)
#inset_axes.plot(time, pegs, 'thistle')
#inset_axes.plot(time,  p_pegs,    'darkred', linestyle='--')
#inset_axes.set_xlim([0, 20])
#inset_axes.set_ylim([-3e-11, -3e-11 + 6e-11])


#axacc.add_patch(Rectangle((0, -3e-11), 20, 6e-11, linewidth=1, edgecolor='k', alpha=0.5, facecolor='none'))



# Multiply pegs by erf taper:
p_pegs_erf = p_pegs * yerf
axacc.plot(time, p_pegs_erf)


# FD of ERF windowed trace:
dx = time[1] - time[0]
acc_2dx = (-p_pegs_erf[4:] + 16 * p_pegs_erf[3:-1] - 30 * p_pegs_erf[2:-2] + 16 * p_pegs_erf[1:-3] - p_pegs_erf[:-4]) / (12 * dx * dx)
ax2nd.plot(time[2:-2],  acc_2dx, 'k')

window = cs(time)

# Plot tapered acceleration trace:
axtap.plot(time, window, 'thistle')
axtap.plot(time, yerf, 'k')



# Plot taper for reference
#axt.plot(time, window, 'darkslategrey')

# Original 2nd derivative with Tukey
axfin.plot(time[2:-2], dd_pegs*window[2:-2], color='thistle', linestyle='-' )

# 4th order quartic without ERF
axfin.plot(time, window*p_pegs_dd, 'darkred', linestyle='--', alpha=0.8)
axfin.plot(time[2:-2], window[2:-2]*acc_2dx, 'k' )


leg_fs = 8

axacc.legend([r'Synthetic acceleration, $\ddot{\Lambda}$ ', 'Quartic polynomial, $\ddot{q}$ '], frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.265, 0.18, 0.4, 0.15) )
axtap.legend(['Tukey Taper', 'ERF taper'], frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.55, 0.48, 0.2, 0.1))
ax2nd.legend(['Synthetic', 'Quartic polynomial', 'Erf-tapered quartic polynomial'], frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.57, 0.3, 0.2, 0.1))
axfin.legend(['Original w. Tukey', 'quartic poly acc w. Tukey', 'ERF-tapered quartic poly acc w.Tukey'], frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.71, 0.3, 0.2, 0.1))
#axt.legend(['Erf taper', 'Tukey taper'])



axacc.tick_params(labelleft=True, labelbottom=False, left=True, bottom=False)
axacc.spines['right'].set_visible(False)
axacc.spines['top'].set_visible(False)
axacc.spines['bottom'].set_visible(False)


#axacc.text(0.95, 0.95, '(a)', transform=axacc.transAxes, fontsize=10, verticalalignment='top')
#inset_axes.text(0.88, 0.95, '(e)', transform=inset_axes.transAxes, fontsize=10, verticalalignment='top')




ax2nd.tick_params(labelleft=True, labelbottom=False, left=True, bottom=False)
ax2nd.spines['top'].set_visible(False)
ax2nd.spines['right'].set_visible(False)
ax2nd.spines['bottom'].set_visible(False)

axtap.tick_params(labelleft=True, labelbottom=False, left=True, bottom=False)
axtap.spines['top'].set_visible(False)
axtap.spines['right'].set_visible(False)
axtap.spines['bottom'].set_visible(False)
axtap.set_ylim([0,1.])


axfin.tick_params(labelleft=True, labelbottom=False, left=True, bottom=False)
axfin.spines['top'].set_visible(False)
axfin.spines['right'].set_visible(False)
axfin.spines['bottom'].set_visible(False)



# Y LABELS:
ylbl_fs = 9

axacc.set_ylabel(r'Acceleration [$m/s^2$]', fontsize=ylbl_fs)
axtap.set_ylabel(r'Taper' + '\nmagnitude', fontsize=ylbl_fs)
ax2nd.set_ylabel(r'2nd deriative of' + '\n  Acceleration [$m/s^4$]', fontsize=ylbl_fs)
axfin.set_ylabel(r'Tapered 2nd deriative' + '\n of Acceleration [$m/s^4$]', fontsize=ylbl_fs)


ytick_fs = 8
v = 2e-10
axacc.set_yticks([-v/2, v])
axacc.set_yticklabels([-v/2, v], fontsize=ytick_fs)

axtap.set_yticks([0, 1])
axtap.set_yticklabels([0, 1], fontsize=ytick_fs)

v = 4e-14
ax2nd.set_yticks([-v, 0, v])
ax2nd.set_yticklabels([-v, 0, v], fontsize=ytick_fs)
ax2nd.set_ylim([-v, v])


axfin.set_yticks([-v, v])
axfin.set_yticklabels([-v, v], fontsize=ytick_fs)
axfin.set_ylim([-v, v])




axxax.tick_params(labelleft=False, left=False)
axxax.spines['top'].set_visible(False)
axxax.spines['right'].set_visible(False)
axxax.spines['left'].set_visible(False)
axxax.set_xlabel('Time [s]', fontsize=ylbl_fs)


tmin = 0
tmax = time[-1]
for axi in [axxax, axacc, ax2nd, axfin, axtap]:
    axi.set_xlim([tmin, tmax])

xlist = np.arange(tmin, tmax, 50).astype(int)
axxax.set_xticks(xlist)
axxax.set_xticklabels(xlist, fontsize=ytick_fs)



fig.set_tight_layout(True)


# Note that we can use the last two of the poly fit at the end of the time series to replace the two that are missing
# by using the 5 point stencil

adj_src  = np.zeros(len(time))
adj_time = np.zeros(len(time))

adj_time[:] = time[:]
adj_src[2:-2]  = window[2:-2]*acc_2dx

# Lets now take the last two values from the polynomial tapered that doesnt have the err function
adj_src[-2:] = window[-2:]*p_pegs_dd[-2:]
adj_src[:2] = 0



for axi in [axxax, axacc, ax2nd, axfin, axtap]:
    axi.axvline(35)


outpath = f'./{adj_path}/adjoint_source_{stn}_{channel}'
np.savetxt(outpath, X=np.array([adj_time, adj_src]))
print(f"saved to {outpath}")

plt.show()
