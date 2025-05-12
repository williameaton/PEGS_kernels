# Script reproduces Figure 2 (MDJ) from Vallee et al 2017: DOI: 10.1126/science.aao0746
import matplotlib.pyplot as plt
import sys
sys.path.append('../classes')
from pegs import create_pegs

# Colours for plotting and station to plot
clrs    = ['#004488', '#DDAA33']
stn       = 'INU'

ddir = '..//data/data_forward_pegs/'

# Load PEGS data using create_pegs
spec = create_pegs(f'{ddir}/PREM_1C/', code='SPECFEM')
spec.load()

# Load real data and cut data to plot
spec.load_real_data(f"{ddir}/Vallee_data/real_data")
spec.cut_to_plot(plusval=35)

# Create figure
fig, ax = plt.subplots(figsize=(8, 5))

# Plot real data
l0, = ax.plot(spec.real['time'][stn], spec.real['pegs'][stn], color='#BB5566', linewidth=3, label=r'Real acc. ($A$)')

# Plot synthetics
l3, = ax.plot(spec.time[stn],  spec.z_grav[stn], clrs[0], linewidth=3, label='Gravity perturb. ($-\mathbf{\delta g}$)')
l3, = ax.plot(spec.time[stn], -spec.z_acc[stn],  clrs[1], linewidth=3, label=r'Ground acc. ($-\partial_t^2 \mathbf{s}$)')
l3, = ax.plot(spec.time[stn],  spec.z_pegs[stn], 'k',     linewidth=3, label='Net PEGS signal ($\partial_t^2 \mathbf{s} - \mathbf{\delta g}$)')

# Set limits
ax.set_xlim([0, 80])
ax.set_ylim([-1.25, 3.5])

# Set y ticks
#ax.set_yticks([-1, 0, 1, 2, 3])
#ax.set_yticklabels(["-1", "0", "1", "2", "3"])

# Set labels
fs = 12
ax.set_xlabel('Time [s]', weight='bold', fontsize=fs)
ax.set_ylabel('Vertical acceleration [nm/s/s]', weight='bold', fontsize=fs)

# Remove ugly spines!
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
legend = ax.legend(fontsize=13, loc='upper left')
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('black')
legend.get_frame().set_alpha(None)

# Event data
ax.text(s='STATION: MDJ', x=0, y=-1.1, fontsize=11, weight='bold')
ax.text(s='M9.0 TŌHOKU-OKI, JAPAN (2011) ', x=40, y=-1.1, fontsize=11, weight='bold')

# Save to pdf and plot to stdout
outname = '/Users/eaton/Documents/Princeton/PEGS_kernels/Figures/pdfs/reproduce_Vallee_INU_example.pdf'
plt.savefig(outname)
plt.show()
print(f'Saved pdf to {outname}')