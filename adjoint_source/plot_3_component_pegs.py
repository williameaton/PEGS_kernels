import obspy
import numpy as np
import matplotlib.pyplot as plt
from wetools import obspy_gen_mpl

fig, ax = plt.subplots(2)


for chl in ['Z', 'N','E']:
    zat, za = obspy_gen_mpl(obspy.read(f'MDJ_DATA/{chl}acc/MDJ_proc.sac')[0])
    zat, zg = obspy_gen_mpl(obspy.read(f'MDJ_DATA/{chl}grav/MDJ_proc.sac')[0])
    ax[0].plot(zat, za+zg)

    t, y = np.loadtxt(f'MDJ_DATA/adjoint_source_MDJ_{chl}')
    ax[1].plot(t,y)


plt.show()