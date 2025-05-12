import matplotlib.pyplot as plt
import sys
sys.path.append('../classes')
from pegs import create_pegs
from wetools.plotting import setup_we_mpl
setup_we_mpl()

# Colours for plotting and station to plot
clrs    = ['#DDAA33' ,'#004488']
stn       = 'MDJ'
fs = 12

ddir = '..//data/data_forward_pegs/'

# Load PEGS data using create_pegs
spec = create_pegs(f'{ddir}/PREM_2C', code='SPECFEM', three_component=True)
spec.load()
spec.cut_to_plot(plusval=35)

# Create figure
fig, ax = plt.subplots(3, figsize=(8, 5), sharex=True)


# Plot synthetics
chls = 'zne'
for ichl in range(3):
    chl = chls[ichl]
    l3, = ax[ichl].plot(spec.time[stn],  spec.traces[f"{chl}_grav"][stn], clrs[0], linewidth=3, label='Gravity perturb. ($\mathbf{\delta g}$)')
    l3, = ax[ichl].plot(spec.time[stn],  spec.traces[f"{chl}_acc"][stn], clrs[1],   linewidth=3, label='Reversed ground acc. ($-\partial_t^2 \mathbf{s}$)')
    l3, = ax[ichl].plot(spec.time[stn],  spec.traces[f"{chl}_pegs"][stn], 'k',     linewidth=3, label='Net PEGS signal ($\partial_t^2 \mathbf{s} - \mathbf{\delta g} $)')

    # Set y ticks
    #ax[ichl].set_yticks([-1.5, -1, -0.5, 0, 0.5])
    #ax[ichl].set_yticklabels(('-1.5', '-1', '-0.5', '0', '0.5'))

    ax[ichl].text(x=5, y = -6.2, s=f'{chl.upper()} channel', weight='bold', fontsize=fs)

    # Remove ugly spines!
    ax[ichl].spines['top'].set_visible(False)
    ax[ichl].spines['right'].set_visible(False)

    # Set labels
    ax[1].set_xlabel('Time [s]', weight='bold', fontsize=fs)


# # Set limits
# ax[0].set_xlim([0, 75])
#
# ax[0].set_ylim([-0.8, 2.2])
# ax[1].set_ylim([-0.3, 2.9])
# ax[2].set_ylim([-0.1, 1.27])


ax[0].set_ylabel('Acceleration [nm/s/s]', weight='bold', fontsize=fs)

# Legend
legend = ax[0].legend(fontsize=9)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)

# Event data
#ax[-1].text(s='STATION: MDJ', x=0, y=-1, fontsize=11, weight='bold')
#ax[-1].text(s='M9.0 TŌHOKU-OKI, JAPAN (2011) ', x=-17, y=0.9, fontsize=11, weight='bold')

# Save to pdf and plot to stdout
outname = '/Users/eaton/Documents/Princeton/PEGS_kernels/forward_simulations/Tohoku/figures/MDJ_allchannels.pdf'
plt.savefig(outname)
plt.show()
print(f'Saved pdf to {outname}')