import obspy
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from wetools import obspy_gen_mpl

# Plots anelastic and elastic to show effect of anelasticity

dir1 = 'elastic_test/anelastic'
dir2 = 'elastic_test/elastic'

stn_list = ['ULN', 'XAN', 'MA2', 'BJT', 'NE93', 'INCN', 'MDJ', 'INU', 'MAJO']

chnls = ['Zacc', 'Grav']

for stn in stn_list:
    fig, ax = plt.subplots(2)

    for i in range(len(chnls)):
        z1 = obspy.read(f"{dir1}/{chnls[i]}/{stn}_proc.sac")
        z2 = obspy.read(f"{dir2}/{chnls[i]}/{stn}_proc.sac")

        t1, z1 = obspy_gen_mpl(z1[0])
        t2, z2 = obspy_gen_mpl(z2[0])

        ax[i].plot(t1, z1)
        ax[i].plot(t2, z2, '--')

        ax[i].set_xlim([0, 350])
plt.show()

