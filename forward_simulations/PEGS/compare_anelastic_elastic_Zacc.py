# Plots anelastic and elastic to show effect of anelasticity
# They are minimal
import obspy
import matplotlib.pyplot as plt
from wetools.funcs import obspy_gen_mpl
from wetools.plotting import setup_we_mpl
import sys
sys.path.append('../../classes')
from pegs import create_pegs

setup_we_mpl()


ddir = '../../data/data_forward_pegs/'
dir1 = f'{ddir}/elastic_test/anelastic'
dir2 = f'{ddir}/elastic_test/elastic'

# Create PEGS object to get the cutoff times:
cutoff = create_pegs(f'{ddir}/qssp/prem_elastic', code='QSSP').cutoff

stn_list = ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO']

chnls = ['Zacc', 'Grav']
for stn in stn_list:
    fig, ax = plt.subplots(2)

    icut = cutoff[stn]

    for i in range(len(chnls)):
        z1 = obspy.read(f"{dir1}/{chnls[i]}/{stn}_proc.sac")
        z2 = obspy.read(f"{dir2}/{chnls[i]}/{stn}_proc.sac")

        t1, z1 = obspy_gen_mpl(z1[0])
        t2, z2 = obspy_gen_mpl(z2[0])

        ax[i].plot(t1[t1 <= icut], z1[t1 <= icut])
        ax[i].plot(t2[t2 <= icut], z2[t2 <= icut], '--')

        ax[i].legend(['Elastic', 'Anelastic'])
        ax[0].set_title(f'STATION: {stn}')

plt.show()

