import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
import scipy.integrate as si
from colour_schemes import hex

# Simulation params:
hdur = 70
M0 = 5.31e22
hdur_gauss = hdur/1.628


# Plotting
highcontrast = hex.Hexes().hex['HighContrast_PT']

# based on specfem timing
t = np.linspace(-115, 115, 3000)


# ------------------------------------ GAUSSIAN STF ------------------------------------

gauss_deriv =   (1/((np.pi**0.5)*hdur_gauss)) * np.exp(- ((t)/hdur_gauss)**2 )
gauss_cum   = si.cumtrapz(y=gauss_deriv, x=t, initial=0)

# ------------------------------------ SIN SQ STF ---------------------------------------
# raw stf
t0 = 70
tmask = np.logical_and(t>=-hdur, t<=hdur)
sin_deriv = t*0
sin_deriv[tmask] =  (1/hdur)*np.sin( np.pi*(t[tmask]-t0)/(2*hdur) )**2
# cumulative
sin_cum = t*0
a = np.pi/(2*hdur)
sin_cum[tmask] =    (1/hdur)*((t[tmask]-t0)/2 - np.sin(2*a*(t[tmask]-t0))/(4*a)  + hdur)
sin_cum[t>=hdur] = 1

# ------------------------------------ ISOSCELES STF ------------------------------------
tri_height = 1/hdur
iso_deriv =  t*0
# Positive gradient
tmask = np.logical_and(t>=-hdur, t<=0)
iso_deriv[tmask] = t[tmask]*(tri_height/hdur) + tri_height
# Negative gradient
tmask = np.logical_and(t<=hdur, t>=0)
iso_deriv[tmask] = t[tmask]*(-tri_height/hdur) + tri_height
iso_cum = si.cumtrapz(y=iso_deriv, x=t, initial=0)




# ------------------------------------ PLOT FIGURE STFs ------------------------------------
fig, ax = plt.subplots(1,2,figsize=(9,5))
fig.set_tight_layout(True)

# Start at 0 instead of -hdur
t+= hdur

# Plot Gaussian, Sinsq, triangle:
mdots     = [gauss_deriv, sin_deriv, iso_deriv]
mdots_cum = [gauss_cum, sin_cum, iso_cum]

# Add curves
for icurve in range(3):
    # raw stf
    ax[0].plot(t, mdots[icurve]/np.amax(mdots[icurve]), highcontrast[icurve+1])
    # cumulative stf
    ax[1].plot(t, mdots_cum[icurve], color=highcontrast[icurve+1])

# ------------------------------------ PLOT BELLS AND WHISTLES -----------------------------
for iax in range(2):
    # Background dash lines
    for ttt in [0, 2*hdur]:
        ax[iax].axvline(ttt, color='k', alpha= 0.1,
                        linestyle='--', zorder=0)

    # Remove spines/labels etc:
    ax[iax].tick_params(labelleft=False, left=False, bottom=False, labelbottom=False)
    ax[iax].spines['top'].set_visible(False)
    ax[iax].spines['right'].set_visible(False)
    ax[iax].spines['left'].set_visible(False)
    ax[iax].spines['bottom'].set_visible(False)

    # y limits
    ax[iax].set_ylim([-0.05, 1.02])

    # Add '2 x half duration' annotation with two arrows pointing opposite directions
    for iar in range(2):
        ax[iax].arrow(iar*2*hdur, -0.04, ((-1)**iar)*2*hdur, 0, head_width=0.015,
                    head_length=3, linewidth=1, color='k', length_includes_head=True)
    ax[iax].text(x=28,y=-0.08,s='2 x half duration', fontsize=11, weight='bold')



# Legend
legend = ax[1].legend(['Gaussian', 'Squared sinusoid', 'Isosceles triangle'], loc=(0.542, 0.05/1.07) )
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)

# Titles
ax[0].set_title('Source time function',      fontsize=11, weight='bold')
ax[1].set_title('Cumulative energy release', fontsize=11, weight='bold')

fig.savefig('./pdfs/compare_stfs.pdf', format='pdf')
plt.show()

