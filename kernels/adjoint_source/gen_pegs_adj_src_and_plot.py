# Plotting the FD time derivatives of the acceleration used in the adjoint source:
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import obspy
from scipy.signal.windows import tukey
from scipy.special import erf
import matplotlib.gridspec as gridspec
from colour_schemes import hex
import numpy as np
from wetools.funcs import obspy_gen_mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def cs(time, alpha):
    # Tukey taper
    window = tukey(len(time), alpha)
    return window

stn       = 'MDJ'
network   = 'IC'
data_fdir =  '../MDJ_kernel/mdj_forward_output_files'

channel   = 'E'
hdur      = 70

# Parameters we choose:
tukey_alpha = 0.15  # Alpha value for Tukey taper
tukey_buffer = 10   # Num of timesteps = 0 buffer on either side of Tukey
poly_deg    = 6     # Polynomial degree for polynomial fit to pegs signal

# Plotting:
ytick_fs = 8        # Y tick font size
ylbl_fs = 9         # Y label font size
polyls   = '--'     # Polynomial trace linestyle
leg_fs   = 8        # Legend font size
txt10_x = 0.02      # Order of magnitude x position
txt10_y = 1         # Order of magnitude y position

# Colours:
HighContrast = hex.Hexes().hex['HighContrast_PT']
pegs_clr  = HighContrast[1]
poly_clr = HighContrast[2]
final_clr = HighContrast[3]


# ---------------------- LOAD PEGS DATA ----------------------
store = [[], []]
types =  [f'{channel}acc', f'{channel}grav']
for i_trtype in range(2):
    trtype = types[i_trtype]
    fpath = f'{data_fdir}/{trtype}/{stn}_proc.sac'

    # Load in the sac:
    trace = obspy.read(fpath)

    # Convert to numpy
    time, acc = obspy_gen_mpl(trace[0])
    dx = np.mean(time[1:] - time[:-1])

    # calculate 5 point finite difference in time
    acc_dx = (-acc[4:] + 8*acc[3:-1] - 8*acc[1:-3] + acc[:-4])/(12*dx)
    acc_2dx =  (-acc[4:] + 16*acc[3:-1] -30*acc[2:-2] + 16*acc[1:-3]
                - acc[:-4])/(12*dx*dx)

    store[i_trtype].append(acc)
    store[i_trtype].append(acc_dx)
    store[i_trtype].append(acc_2dx)

# Combine ground acceleration and gravitational acceleration:
pegs    = store[1][0] + store[0][0]
d_pegs  = store[1][1] + store[0][1]
dd_pegs = store[1][2] + store[0][2]

# Timestep data:
ntimesteps = len(time)
dx = time[1] - time[0]

# ---------------------- POLYNOMIAL FIT TO PEGS DATA ----------------------
p = np.polyfit(x=time, y=pegs, deg=poly_deg)

# Compute polynomial pegs signal and its derivatives
p_pegs    = np.zeros(ntimesteps)
p_pegs_dd = np.zeros(ntimesteps)
for id in range(poly_deg+1):
    p_pegs +=  p[id]*(time**(poly_deg-id))
for id in range(poly_deg-1):
    p_pegs_dd += (poly_deg-id)*(poly_deg-1-id)*p[id]*(time**(poly_deg-2-id))

# ----------------------TAPER POLYNOMIAL SIGNAL ----------------------

# Create error function
xerf = np.linspace(-4.5, 12, len(time))
yerf = 0.5 + 0.5 * erf(xerf)

# Multiply polynomiall pegs by erf taper:
p_pegs_erf = p_pegs * yerf

# FD of ERF-windowed trace:
erf_poly_2nd = (-p_pegs_erf[4:] + 16 * p_pegs_erf[3:-1] - 30 * p_pegs_erf[2:-2] +
                16 * p_pegs_erf[1:-3] - p_pegs_erf[:-4]) / (12 * dx * dx)

# Create Tukey taper
window = np.zeros(ntimesteps)
window[tukey_buffer:-tukey_buffer] = cs(time[tukey_buffer:-tukey_buffer], tukey_alpha)



# ---------------------- SAVE DATA FOR USE IN ADJ SIMULATION ----------------------
# Note that we can use the last two of the poly fit at the end of the time series
# to replace the two that are missing by using the 5 point stencil
adj_src  = np.zeros(ntimesteps)
adj_time = np.zeros(ntimesteps)

adj_time[:] = time[:] - hdur
adj_src[2:-2]  = window[2:-2]*erf_poly_2nd

# Lets now take the last two values from the polynomial tapered that doesnt have the err function
adj_src[-2:] = 0
adj_src[:2] = 0
np.savetxt(f'{data_fdir}/adj_source/{stn}.{network}.MX{channel}.adj', X=np.array([adj_time, adj_src]).T)


# ---------------------- PLOT FIGURE ----------------------
fig = plt.figure(figsize=(4,7.3))
fig.set_tight_layout(True)

spec = gridspec.GridSpec(ncols=1, nrows=15, figure=fig)
axacc = fig.add_subplot(spec[0:4, 0])       # Acceleration
axtap = fig.add_subplot(spec[4:6, 0])       # Taper
ax2nd = fig.add_subplot(spec[6:10, 0])      # 2nd derivative of acceleration
axfin = fig.add_subplot(spec[10:14, 0])     # Final tapered result
axxax = fig.add_subplot(spec[-1, 0])        # Separated x axis

axeslist = [axacc, axtap, ax2nd, axfin, axxax]     # ignoring the x axis (axxax)


# Plot SPECFEM pegs signal
axacc.plot(time, pegs,   color=pegs_clr)                    # original pegs acceleration from SPECFEM
axacc.plot(time, p_pegs, color=poly_clr, linestyle=polyls)  # polynomial fit to pegs acceleration

# Plot 2nd derivative traces
ax2nd.plot(time[2:-2], dd_pegs, color=pegs_clr)             # 2nd derivative of SPECFEM pegs signal
ax2nd.plot(time,  p_pegs_dd, poly_clr, linestyle=polyls)    # polynomial trace 2nd derivative without ERF
ax2nd.plot(time[2:-2],  erf_poly_2nd, color=final_clr)      # plot Erf-windowed polynomial trace 2nd derivative

# Plot inset zoom in of the quartic vs PEGS acceleration in acceleration subplot
inset_axes = inset_axes(axacc,
                        loc = 'center left',
                        width="50%",
                        height = "50%",
                        bbox_to_anchor=(0.03, 0.4, 0.8, 0.5),
                        bbox_transform=axacc.transAxes
                        )
inset_axes.tick_params(labelleft=False, labelbottom=False, left=False, bottom=False)
inset_axes.plot(time, pegs, pegs_clr)
inset_axes.plot(time,  p_pegs,    poly_clr, linestyle=polyls)
inset_axes.set_xlim([0, 20])
inset_axes.set_ylim([-0.2e-11, 0.2e-11])
for axi in ['right','left','top', 'bottom']:
    inset_axes.spines[axi].set_alpha(0.5)
axacc.add_patch(Rectangle((0, -3e-11), 20, 6e-11, linewidth=1, edgecolor='k', alpha=0.5, facecolor='none'))

# Plot tapered acceleration trace:
axtap.plot(time, window, color=pegs_clr)    # Tukey window taper
axtap.plot(time, yerf,   color=final_clr)   # Erf taper

# Plot 2nd derivatives with window (Tukey) taper added
axfin.plot(time[2:-2], dd_pegs*window[2:-2], color=pegs_clr)                # original pegs signal
axfin.plot(time, window*p_pegs_dd, poly_clr, linestyle=polyls, alpha=0.8)   # polynomial + tukey taper
axfin.plot(time[2:-2], window[2:-2]*erf_poly_2nd, color=final_clr )         # polynomial + tukey + erf taper

# Legends
axacc.legend([r'Synthetic acceleration, $\ddot{\Lambda}$ ', 'Quartic polynomial, $\ddot{q}$ '],
             frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.225, 0.12, 0.4, 0.15) )
axtap.legend(['Tukey taper', 'ERF taper'],
             frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.75, 0.45, 0.2, 0.1))
ax2nd.legend(['Synthetic', 'Quartic polynomial', 'Erf-tapered quartic polynomial'],
             frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.52, 0.245, 0.2, 0.1))
axfin.legend(['Original with Tukey', 'Quartic polynomial\nacceleration with Tukey',
              'ERF-tapered quartic polynomial\nacceleration with Tukey'],
             frameon=False,  fontsize=leg_fs,  bbox_to_anchor=(0.545, 0.41, 0.2, 0.1))


YTICKS     = [[-1e-9, 0],        [0, 1],       [-2e-12, 0, 2e-12], [-1e-12, 0]      ]
TICKLABELS = [[-1, '0'],         [0, 1],       [-2, 0, 2],         [-1, 0]          ]
YLIMS      = [[-1.4e-9, 0.2e-9], [-0.02,1.02], [-3e-12, 3e-12],    [-1e-12, 0.2e-12]]
LABELS     = [r'Acceleration [$\mathbf{ms}^\mathbf{-2}$]',
              r'Taper' + '\nmagnitude',
              r'2$^{\mathbf{nd}}$ derivative of' + '\n  acceleration [$\mathbf{m s}^\mathbf{-4}$]',
              r'Tapered 2$^{\mathbf{nd}}$ derivative' + '\n of acceleration [$\mathbf{m s}^\mathbf{-4}$]'
             ]


iax = 0
for axi in axeslist[:-1]:
    # Remove ticks and spines
    axi.tick_params(labelleft=True, labelbottom=False, left=True, bottom=False)
    for spine in ['right', 'top', 'bottom']:
        axi.spines[spine].set_visible(False)
    # Set tick locations, strings, ylims and labels
    axi.set_yticks(YTICKS[iax])
    axi.set_yticklabels(TICKLABELS[iax], fontsize=ytick_fs)
    axi.set_ylim(YLIMS[iax])
    axi.set_ylabel(LABELS[iax], fontsize=ylbl_fs, weight='bold')
    iax += 1


# Add text showing the order of magnitude for the y axes
yax_order10 = [9, 12, 12]
iax = 0
for axi in [axacc, ax2nd, axfin]:
    axi.text(txt10_x, txt10_y, 'x $10^{-' + str(yax_order10[iax]) + '}$' ,  transform=axi.transAxes, fontsize=8, verticalalignment='top')
    iax +=1

# Align the y axes
for iax in axeslist[:-1]:
    iax.yaxis.set_label_coords(-0.1, 0.5)

# Set the x limits
for axi in axeslist:
    axi.set_xlim([0, time[-1]])

# Add x axis that is offset from other plots (personal preference)
axxax.tick_params(labelleft=False, left=False)
for spine in ['top', 'right', 'left']:
    axxax.spines[spine].set_visible(False)
axxax.set_xlabel('Time [s]', fontsize=ylbl_fs, weight='bold')
axxax.set_xticks([0, 50, 100, 150])
axxax.set_xticklabels([0, 50, 100, 150], fontsize=ytick_fs)


fig.savefig(f"../../Figures/pdfs/example_adjoint_source_pegs.pdf", format='pdf')
plt.show()

